source("classes/utils.r")

library(dplyr)
library(cluster)
# library(umap) removed: nothing in this class calls it. RunUMAP below is
# Seurat's, which uses uwot. The R `umap` package is not installed in the
# scRNA_seq env and this import killed a job at source() time.
library(pheatmap)
library(Hmisc)
library(Seurat)
library(plyr)
library(ggplot2)
library(scales)

#' Clusterer Class for Clustering Seurat Objects
#'
#' This R6 class provides methods to perform clustering on Seurat objects,
#' identify the best clustering resolution using silhouette scores, set the
#' best clusters, and find top marker genes for each cluster.
#'
#' @field seurat_obj A Seurat object containing the integrated single-cell
#'        data.
#' @field resolutions A vector of clustering resolutions to evaluate.
#' @field dims A vector of PCA dimensions to use for clustering.
#' @field verbose A logical value indicating whether to print verbose output.
Clusterer <- R6Class( # nolint
  "Clusterer",
  public = list(
    seurat_obj = NULL,
    verbose = NULL,
    resolutions = NULL,
    dims = NULL,
    seed = 1234,
    utils = Utils$new(),

    #' Initialize the Clusterer
    #'
    #' @param seurat_obj  A Seurat object containing the integrated single-cell
    #'                    data.
    #' @param verbose     A logical value indicating whether to print verbose
    #'                    output (default: FALSE).
    #' @param resolutions A vector of clustering resolutions to evaluate
    #'                    (default: seq(0.1, 1, by = 0.1)).
    #' @param dims        A vector of PCA dimensions to use for clustering
    #'                    (default: 1:50).
    #' @param seed        RNG seed, set before the resolution sweep and before
    #'                    the silhouette subsampling.
    #'
    #'                    WITHOUT THIS THE CLUSTER COUNT IS NOT REPRODUCIBLE.
    #'                    calc_silhouette_scores draws 100 subsamples of 1,000
    #'                    cells PER RESOLUTION to score each one, via an
    #'                    unseeded sample(). best_resolution is then the argmax
    #'                    over those noisy estimates, so two runs of identical
    #'                    code on identical input could select different
    #'                    resolutions and therefore report a different number of
    #'                    clusters. Defaults to 1234.
    initialize = function(seurat_obj, verbose = FALSE,
                          resolutions = seq(0.1, 1, by = 0.1), dims = 1:50,
                          seed = 1234) {
      self$seurat_obj <- seurat_obj
      self$verbose <- verbose
      self$resolutions <- resolutions
      self$dims <- dims
      self$seed <- seed
    },

    #' Run Clustering
    #'
    #' This method performs PCA, UMAP, neighbor finding, and clustering on the
    #' Seurat object.
    #'
    #' @return The Seurat object with clustering results.
    run_clustering = function() {
      set.seed(self$seed)
      self$seurat_obj <-
        RunPCA(self$seurat_obj,
               features = VariableFeatures(object = self$seurat_obj),
               verbose = self$verbose)
      self$seurat_obj <-
        RunUMAP(self$seurat_obj, dims = self$dims, verbose = self$verbose)
      self$seurat_obj <-
        FindNeighbors(self$seurat_obj, dims = self$dims,
                      verbose = self$verbose)
      self$seurat_obj <-
        FindClusters(self$seurat_obj, resolution = self$resolutions,
                     verbose = self$verbose, algorithm = 1)
      return(self$seurat_obj)
    },

    #' Find the Best Clustering Resolution Using Silhouette Scores
    #'
    #' This function identifies the best clustering resolution for a Seurat
    #' object by calculating silhouette scores for various resolutions and
    #' selecting the one with the highest mean silhouette score.
    #'
    #' @param resolutions A vector of clustering resolutions to evaluate.
    #' @param pca_dims A vector of PCA dimensions to use for silhouette score
    #'                 calculation.
    #'
    #' @return A list containing:
    #'  - best_resolution: The resolution with the highest mean silhouette
    #'                     score.
    #'  - mean_scores:     A named vector of mean silhouette scores for each
    #'                     resolution.
    #'  - sd_scores:       A named vector of standard deviations of silhouette
    #'                     scores for each resolution.
    calc_silhouette_scores = function() {
      set.seed(self$seed)
      # FindClusters names its columns "<DefaultAssay>_snn_res.<r>", so this
      # prefix MUST follow the assay. Hardcoding "integrated_snn_res." is right
      # in Step 4 and silently wrong in Step 7, where DefaultAssay is "VIPER":
      # FindClusters writes VIPER_snn_res.*, while this read the STALE Step-4
      # columns. Job 60551 therefore returned the GENE-EXPRESSION clustering as
      # its VIPER result - same 19 clusters, same sizes, same Origin
      # composition - and scored silhouettes for gene clusters against VIPER
      # PCA coordinates. Job 60626 did the same on the macrophage arm.
      assay_prefix <- paste0(DefaultAssay(self$seurat_obj), "_snn_res.")
      clustering_results <-
        self$seurat_obj@meta.data[, grepl(assay_prefix,
                                          colnames(self$seurat_obj@meta.data),
                                          fixed = TRUE), drop = FALSE]
      if (ncol(clustering_results) == 0L)
        stop("No clustering columns found matching '", assay_prefix, "'")
      pca_matrix <-
        as.data.frame(
          # [, self$dims] selects PCs. This was [self$dims, ], which selects
          # the first 50 CELLS instead - so every silhouette sweep before this
          # fix chose best_resolution from 50 cells rather than 50 PCs.
          t(Embeddings(self$seurat_obj, reduction = "pca")[, self$dims])
        )
      silhouette_scores <-
        private$compute_silhouette_scores(pca_matrix, clustering_results)
      mean_scores <- silhouette_scores$means
      sd_scores <- silhouette_scores$sd

      best_resolution <-
        tail(self$resolutions[which(mean_scores == max(mean_scores))], n = 1)

      return(list(best_resolution = best_resolution, mean_scores = mean_scores,
                  sd_scores = sd_scores))
    },

    #' Set Clusters Based on the Best Resolution
    #'
    #' This function assigns clusters to a Seurat object based on the best
    #' clustering resolution and sets the active identity class to these
    #' clusters.
    #'
    #' @param best_resolution The best clustering resolution identified.
    #'
    #' @return                The Seurat object with clusters set based on the
    #'                        best resolution.
    set_best_clusters = function(best_resolution) {
      cluster_column_name <-
        paste0(DefaultAssay(self$seurat_obj), "_snn_res.", best_resolution)
      if (!cluster_column_name %in% colnames(self$seurat_obj@meta.data)) {
        stop("No column '", cluster_column_name, "' in meta.data. Available: ",
             paste(grep("_snn_res[.]", colnames(self$seurat_obj@meta.data),
                        value = TRUE), collapse = ", "))
      }
      self$seurat_obj$seurat_clusters <-
        self$seurat_obj@meta.data[, cluster_column_name]
      Idents(self$seurat_obj) <- "seurat_clusters"
      return(self$seurat_obj)
    },

    #' Find top genes for each cluster
    #'
    #' This function identifies the top marker genes for each cluster in the
    #' integrated Seurat object. It uses the scaled gene expression data to
    #' perform differential expression analysis between each cluster and all
    #' other cells, and selects the top genes based on the specified log fold
    #' change threshold.
    #'
    #' @param assay_name      The name of the assay to use for finding markers.
    #' @param n_top_genes     The number of top genes to find for each cluster.
    #' @param logfc_threshold The log fold change threshold for marker genes.
    #'
    #' @return                A data frame containing the top marker genes for
    #'                        each cluster. The data frame includes the
    #'                        following columns:
    #'                        - `gene`:       The gene name.
    #'                        - `cluster`:    The cluster for which the gene is
    #'                                        a marker.
    #'                        - `avg_log2FC`: The average log2 fold change of
    #'                                        the gene in the cluster compared
    #'                                        to all other cells.
    #' Find Top Master Regulators per Cluster (VIPER / NES data)
    #'
    #' The NES counterpart of find_top_genes, following
    #' scripts/archive/viper_global.R rather than this class's homemade
    #' difference-of-means.
    #'
    #' Three parameters are load-bearing and none of them are Seurat defaults:
    #'   mean.fxn = rowMeans  Seurat's default applies expm1(), which assumes
    #'                        log1p-transformed expression. NES is negative for
    #'                        roughly half of all entries, so expm1 on it is
    #'                        meaningless. Overriding to a plain mean gives the
    #'                        difference in mean NES.
    #'   fc.name  = avg_diff  because the quantity is a difference in NES, not a
    #'                        log fold change. Naming it avg_log2FC invites it to
    #'                        be read and thresholded as one.
    #'   min.pct = 0,         "percent expressed" and a fold-change floor are
    #'   logfc.threshold = 0  expression concepts. A regulator is not "expressed"
    #'                        in a cell; it has an activity. Any nonzero cut here
    #'                        silently drops regulators that differ modestly but
    #'                        consistently.
    #'
    #' @param n_top_genes         Regulators to keep per cluster. Defaults to 5.
    #' @param max_cells_per_ident Cells sampled per cluster for the test.
    #'                            Defaults to 500, as in the reference.
    #' @param assay_name          Defaults to DefaultAssay.
    #'
    #' @return A list with `top` (n_top_genes per cluster, deduplicated) and
    #'         `all` (the full marker table, worth writing to disk).
    find_top_regulators = function(n_top_genes = 5, max_cells_per_ident = 500,
                                   assay_name = NULL) {
      if (is.null(assay_name)) assay_name <- DefaultAssay(self$seurat_obj)
      markers <- FindAllMarkers(self$seurat_obj,
                                assay = assay_name,
                                slot = "data",
                                only.pos = TRUE,
                                min.pct = 0,
                                logfc.threshold = 0,
                                test.use = "t",
                                max.cells.per.ident = max_cells_per_ident,
                                random.seed = 1234,
                                mean.fxn = function(x) rowMeans(x),
                                fc.name = "avg_diff",
                                verbose = self$verbose)
      if (nrow(markers) == 0 || !"cluster" %in% colnames(markers)) {
        stop("FindAllMarkers returned no markers - with ",
             nlevels(Idents(self$seurat_obj)),
             " ident(s) there is nothing to contrast.")
      }
      top <- markers %>%
        dplyr::group_by(cluster) %>%
        dplyr::slice_max(avg_diff, n = n_top_genes, with_ties = FALSE) %>%
        dplyr::distinct(gene, .keep_all = TRUE) %>%
        dplyr::ungroup()
      message(sprintf("  markers: %d regulators over %d cluster(s); keeping %d",
                      nrow(markers), nlevels(factor(markers$cluster)),
                      nrow(top)))
      return(list(top = top, all = markers))
    },

    find_top_genes = function(assay_name = NULL, n_top_genes = 5,
                              logfc_threshold = 0.25) {
      if (is.null(assay_name)) assay_name <- DefaultAssay(self$seurat_obj)
      scale_data <-
        GetAssayData(self$seurat_obj, assay = assay_name, layer = "scale.data")
      clusters <- Idents(self$seurat_obj)

      all_markers <- lapply(unique(clusters), function(cluster) {
        private$find_cluster_markers(scale_data, clusters, cluster,
                                     logfc_threshold, n_top_genes)
      })

      all_markers_df <- bind_rows(all_markers)

      return(all_markers_df)
    }
  ),

  #############################################################################
  #                           PRIVATE METHODS                                 #
  #############################################################################
  private = list(
    #' Compute silhouette scores for multiple Louvain clusterings
    #'
    #' This function evaluates multiple Louvain clusterings of a single-cell
    #' matrix by computing the silhouette scores for a subsample of the data.
    #' It iterates over 100 alternative resolution values, each time
    #' sub-sampling 1000 cells (or fewer if fewer are available) and calculates
    #' the mean and standard deviation of the silhouette scores for each
    #' clustering resolution.
    #'
    #' @param mat   Matrix with rows as principal component vectors and columns
    #'              as samples.
    #' @param clust Matrix with rows as samples and each column as a clustering
    #'              vector for a given resolution.
    #'
    #' @return      List containing the means and standard deviations of
    #'              silhouette scores for each clustering resolution.
    compute_silhouette_scores = function(mat, clust) {
      num_resolutions <- ncol(clust)
      num_subsamples <- 100

      silhouette_scores <-
        private$initialize_silhouette_scores(num_subsamples,
                                             num_resolutions)

      for (resolution_index in 1:num_resolutions) {
        for (subsample_index in 1:num_subsamples) {
          sampled_indices <- private$sample_cells(mat, 1000)
          distance_matrix <-
            self$utils$compute_distance_matrix(mat[, sampled_indices])
          silhouette_scores[subsample_index, resolution_index] <-
            private$compute_silhouette_width(
              clust[sampled_indices, resolution_index],
              distance_matrix
            )
        }
      }

      list(means = colMeans(silhouette_scores, na.rm = TRUE),
           sd = apply(silhouette_scores, 2, sd, na.rm = TRUE))
    },

    #' Initialize matrix for storing silhouette scores
    #'
    #' @param num_subsamples  Number of subsamples to compute.
    #' @param num_resolutions Number of resolution values/clustering vectors.
    #'
    #' @return                Initialized matrix for storing silhouette scores.
    initialize_silhouette_scores = function(num_subsamples, num_resolutions) {
      matrix(rep(NA, num_subsamples * num_resolutions), nrow = num_subsamples)
    },

    #' Randomly sample cells from the matrix
    #'
    #' @param mat       Data matrix.
    #' @param num_cells Number of cells to sample.
    #'
    #' @return          Indices of sampled cells.
    sample_cells = function(mat, num_cells) {
      sample(seq_len(ncol(mat)), min(num_cells, ncol(mat)))
    },

    #' Compute the mean silhouette width for a clustering
    #'
    #' @param clustering      Clustering vector for the subsampled cells.
    #' @param distance_matrix Distance matrix for the subsampled cells.
    #'
    #' @return                Mean silhouette width for the clustering.
    compute_silhouette_width = function(clustering, distance_matrix) {
      if (length(unique(clustering)) <= 1) return(0)

      silhouette_scores <- silhouette(as.numeric(clustering), distance_matrix)
      mean(silhouette_scores[, "sil_width"])
    },

    #' Find Cluster Markers
    #'
    #' This function identifies the top marker genes for a given cluster using
    #' the log fold change (logFC) values calculated from the scaled gene
    #' expression data.
    #'
    #' @param scale_data      A matrix of scaled gene expression data, where
    #'                        rows are genes and columns are cells.
    #' @param clusters        A factor or vector of cluster identities for each
    #'                        cell.
    #' @param cluster         The cluster identity for which to find marker
    #'                        genes.
    #' @param logfc_threshold The log fold change threshold for selecting
    #'                        marker genes.
    #' @param n_top_genes     The number of top marker genes to select for the
    #'                        cluster.
    #'
    #' @return                A data frame containing the top marker genes for
    #'                        the specified cluster.
    find_cluster_markers = function(scale_data, clusters, cluster,
                                    logfc_threshold, n_top_genes) {
      cluster_cells <- clusters == cluster
      other_cells <- clusters != cluster
      logfc <- private$calculate_logfc(scale_data, cluster_cells, other_cells)

      markers <- data.frame(
        gene = rownames(scale_data),
        cluster = cluster,
        avg_log2FC = logfc,
        stringsAsFactors = FALSE
      )

      markers <- markers %>% filter(avg_log2FC > logfc_threshold)

      top_markers <- markers %>% top_n(n = n_top_genes, wt = avg_log2FC)

      return(top_markers)
    },

    #' Calculate Log Fold Change
    #'
    #' This function calculates the average expression and log fold change
    #' (logFC) for each gene between cells in a given cluster and all other
    #' cells.
    #'
    #' @param scale_data    A matrix of scaled gene expression data, where rows
    #'                      are genes and columns are cells.
    #' @param cluster_cells A logical vector indicating which cells belong to
    #'                      the current cluster.
    #' @param other_cells   A logical vector indicating which cells belong to
    #'                      all other clusters.
    #'
    #' @return              A numeric vector of log fold changes for each gene.
    calculate_logfc = function(scale_data, cluster_cells, other_cells) {
      avg_exp_cluster <- rowMeans(scale_data[, cluster_cells, drop = FALSE])
      avg_exp_other <- rowMeans(scale_data[, other_cells, drop = FALSE])
      logfc <- avg_exp_cluster - avg_exp_other
      return(logfc)
    }
  )
)
