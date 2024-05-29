library(dplyr)
library(cluster)
library(umap)
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
#' @field stepping A vector of clustering resolutions to evaluate.
#' @field dims A vector of PCA dimensions to use for clustering.
#' @field verbose A logical value indicating whether to print verbose output.
Clusterer <- R6Class( # nolint
  "Clusterer",
  public = list(
    seurat_obj = NULL,
    stepping = NULL,
    dims = NULL,
    verbose = NULL,

    #' Initialize the Clusterer
    #'
    #' @param seurat_obj A Seurat object containing the integrated single-cell
    #'                   data.
    #' @param stepping   A vector of clustering resolutions to evaluate
    #'                   (default: seq(0.1, 1, by = 0.1)).
    #' @param dims       A vector of PCA dimensions to use for clustering
    #'                   (default: 1:50).
    #' @param verbose    A logical value indicating whether to print verbose
    #'                   output (default: FALSE).
    initialize = function(seurat_obj, stepping = seq(0.1, 1, by = 0.1),
                          dims = 1:50, verbose = FALSE) {
      self$seurat_obj <- seurat_obj
      self$stepping <- stepping
      self$dims <- dims
      self$verbose <- verbose
    },

    #' Run Clustering
    #'
    #' This method performs PCA, UMAP, neighbor finding, and clustering on the
    #' Seurat object.
    #'
    #' @return The Seurat object with clustering results.
    run_clustering = function() {
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
        FindClusters(self$seurat_obj, resolution = self$stepping,
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
      clustering_results <-
        self$seurat_obj@meta.data[, grepl("integrated_snn_res.",
                                          colnames(self$seurat_obj@meta.data))]
      pca_matrix <-
        as.data.frame(
          t(self$seurat_obj@reductions$pca@cell.embeddings[self$dims, ])
        )
      silhouette_scores <-
        private$compute_silhouette_scores(pca_matrix, clustering_results)
      mean_scores <- silhouette_scores$means
      sd_scores <- silhouette_scores$sd

      best_resolution <-
        tail(self$stepping[which(mean_scores == max(mean_scores))], n = 1)

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
        paste("integrated_snn_res.", best_resolution, sep = "")
      self$seurat_obj$seurat_clusters <-
        self$seurat_obj@meta.data[, cluster_column_name]
      Idents(self$seurat_obj) <- "seurat_clusters"
      return(self$seurat_obj)
    },

    #' Find top genes for each cluster
    #'
    #' This function identifies the top marker genes for each cluster in the
    #' integrated Seurat object.
    #'
    #' @param assay_name      The name of the assay to use for finding markers.
    #' @param n_top_genes     The number of top genes to find for each cluster.
    #' @param logfc_threshold The log fold change threshold for marker genes.
    #'
    #' @return                A data frame containing the top marker genes for
    #'                        each cluster.
    find_top_genes = function(assay_name = "SCT", n_top_genes = 10,
                              logfc_threshold = 0.25) {
      # Prepare the SCT assay for differential expression analysis
      self$seurat_obj <- PrepSCTFindMarkers(self$seurat_obj)

      # Find all markers
      all_markers <- FindAllMarkers(
        self$seurat_obj,
        assay = assay_name,
        only.pos = TRUE,
        logfc.threshold = logfc_threshold
      )

      # Check if the cluster column exists
      if (!"cluster" %in% colnames(all_markers)) {
        stop("The 'cluster' column is missing in the markers data frame.")
      }

      # Select top markers for each cluster
      top_genes <- all_markers %>%
        group_by(cluster) %>% # nolint
        top_n(n = n_top_genes, wt = avg_log2FC) # nolint

      return(top_genes)
    }
  ),

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
            private$compute_distance_matrix(mat[, sampled_indices])
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

    #' Compute Distance Matrix
    #'
    #' Calculates a distance matrix using Pearson correlation.
    #'
    #' @param dat_mat A matrix of gene expression data (genes x samples).
    #'
    #' @return        A distance matrix.
    compute_distance_matrix = function(dat_mat) {
      if (!is.matrix(dat_mat)) {
        dat_mat <- as.matrix(dat_mat)
      }

      dist_mat <- as.dist(1 - cor(dat_mat, method = "pearson"))
      return(dist_mat)
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
    }
  )
)
