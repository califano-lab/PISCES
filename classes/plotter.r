library(R6)
library(Seurat)
library(ggplot2)
library(plyr)

#' Plotter Class for Visualizing Clustering Results
#'
#' This class provides methods to plot silhouette scores, UMAP clusters,
#' UMAP with refined labels, and gene heatmaps. It utilizes a Seurat object
#' containing clustering results and allows customization of clustering labels
#' and plotting parameters.
#'
#' @field seurat_obj A Seurat object containing clustering and UMAP results.
#' @field cluster_labels A vector of cluster labels to assign to clusters.
#' @field plot_output_path The directory path where the plots will be saved.
#' @field resolutions A vector of clustering resolutions to evaluate.
Plotter <- R6Class( # nolint
  "Plotter",
  public = list(
    seurat_obj = NULL,
    cluster_labels = NULL,
    plot_output_path = NULL,
    resolutions = NULL,

    #' Initialize the Plotter Object
    #'
    #' @param seurat_obj       A Seurat object containing clustering and UMAP
    #'                         results.
    #' @param plot_output_path The directory path where the plots will be saved.
    #' @param cluster_labels   A vector of cluster labels to assign to clusters.
    #'                         Defaults to NULL.
    #' @param resolutions      A vector of clustering resolutions to evaluate.
    #'                         Defaults to seq(0.1, 1, by = 0.1).
    initialize = function(seurat_obj, plot_output_path, cluster_labels = NULL,
                          resolutions = seq(0.1, 1, by = 0.1)) {
      self$seurat_obj <- seurat_obj
      self$cluster_labels <- cluster_labels
      self$plot_output_path <- plot_output_path
      self$resolutions <- resolutions
    },

    #' Plot Violin Plot
    #'
    #' This function plots a violin plot for a given feature and groups the
    #' data by a specified metadata column. It saves the plot to the specified
    #' output path.
    #'
    #' @param features A character vector of features to plot.
    #' @param group_by A metadata column to group the data by.
    #' @param pt_size  The size of the points in the plot. Defaults to 0.
    #'
    #' @return         None. The function saves the plot as a PNG file in the
    #'                 specified directory.
    plot_violin = function(features, group_by, pt_size = 0) {
      p <-
        VlnPlot(self$seurat_obj, features, group.by = group_by,
                pt.size = pt_size)
      ggsave(file.path(self$plot_output_path,
                       paste0(features, "_violin_plot.png")), plot = p)
    },

    #' Plot Silhouette Scores
    #'
    #' This function plots the mean silhouette scores with error bars for
    #' different clustering resolutions and highlights the best resolution
    #' based on the highest mean silhouette score.
    #'
    #' @param mean_scores A vector of mean silhouette scores for each
    #'                    resolution.
    #' @param sd_scores   A vector of standard deviations of silhouette scores
    #'                    for each resolution.
    #'
    #' @return            None. The function saves the plot as a PNG file in
    #'                    the specified directory.
    plot_silhouette_scores = function(mean_scores, sd_scores) {
      # Open the PNG device FIRST and draw into it, rather than drawing to the
      # default device and dev.copy()ing afterwards.
      #
      # The old order silently produced nothing under Rscript. Non-interactive R
      # has no screen device, so errbar() opened the fallback pdf and every
      # silhouette plot this project ever made accumulated in a stray Rplots.pdf
      # in the working directory - silhouette_scores.png existed nowhere in the
      # output tree. The per-resolution scores were therefore unrecoverable after
      # a run, since resume_from_step6.r logs only the winning resolution, and
      # each configuration had to be re-scored from its saved object just to see
      # its own sweep.
      #
      # on.exit(add = TRUE) guarantees the device closes even if errbar or legend
      # errors, which would otherwise leave a truncated PNG and a dangling device.
      plot_path <- file.path(self$plot_output_path, "silhouette_scores.png")
      png(plot_path, width = 1000, height = 750, res = 120)
      on.exit(dev.off(), add = TRUE)

      errbar(self$resolutions, mean_scores, mean_scores + sd_scores,
             mean_scores - sd_scores, ylab = "Mean Silhouette Score",
             xlab = "Resolution Parameter")
      lines(self$resolutions, mean_scores)

      # na.rm: a resolution that collapsed to a single cluster scores NA, not 0.
      best_resolution <-
        tail(self$resolutions[which(mean_scores == max(mean_scores, na.rm = TRUE))], n = 1)
      legend("topright",
             paste("Best Resolution", best_resolution, sep = " = "))

      # The numbers themselves, so a run records its own sweep and no re-score is
      # needed to read it back.
      message("  silhouette sweep -> ", plot_path)
      for (i in seq_along(self$resolutions))
        message(sprintf("    res %-4s  %.4f +/- %.4f", self$resolutions[i],
                        mean_scores[i], sd_scores[i]))
      invisible(plot_path)
    },

    #' Plot UMAP Clusters
    #'
    #' This function plots UMAP clusters for a Seurat object, assigning
    #' provided cluster labels and saving the plot to the specified output
    #' path.
    #'
    #' @return The Seurat object with updated cluster labels.
    plot_umap_clusters = function() {
      # levels(droplevels(...)), NOT unique(). unique() on a factor returns
      # values in ORDER OF FIRST APPEARANCE, so "Cluster1..N" was assigned by
      # whichever cluster happened to own the first cell - making the labels
      # non-deterministic across runs and impossible to join reliably against
      # the numeric seurat_clusters used in viper_cluster_markers.csv. Since
      # plot_umap_clusters WRITES these labels back over seurat_clusters
      # (mapvalues below), that mislabels the saved object, not just the figure.
      unique_clusters <- levels(droplevels(factor(self$seurat_obj$seurat_clusters)))
      num_clusters <- length(unique_clusters)

      # Ensure the number of provided labels matches the number of clusters
      if (num_clusters > length(self$cluster_labels)) {
        warning(paste0("More clusters found in data than provided labels. ",
                       "Some clusters will be unlabeled."))
        self$cluster_labels <-
          c(self$cluster_labels,
            paste0("Cluster", (length(self$cluster_labels) + 1):num_clusters))
      } else if (num_clusters < length(self$cluster_labels)) {
        self$cluster_labels <- self$cluster_labels[1:num_clusters]
      }

      label_mapping <- setNames(self$cluster_labels, unique_clusters)
      self$seurat_obj@meta.data$seurat_clusters <-
        plyr::mapvalues(self$seurat_obj@meta.data$seurat_clusters,
                        from = unique_clusters, to = self$cluster_labels)

      p <-
        DimPlot(self$seurat_obj, reduction = "umap",
                group.by = "seurat_clusters", label = TRUE, label.size = 7,
                repel = TRUE) + NoLegend()
      umap_cluster_plot_path <-
        file.path(self$plot_output_path, "umap_clustering_results.png")
      ggsave(umap_cluster_plot_path, plot = p, width = 10, height = 8)

      return(self$seurat_obj)
    },

    #' Plot UMAP with Refined Labels
    #'
    #' This function plots a UMAP visualization of the Seurat object with
    #' refined labels.
    #'
    #' @return The UMAP plot is saved to the specified directory.
    plot_umap_with_labels = function() {
      private$filter_blueprint_labels()
      p <- DimPlot(self$seurat_obj, reduction = "umap", label = TRUE,
                   repel = TRUE, label.size = 5,
                   group.by = "refined_labels") + NoLegend()
      ggsave(file.path(self$plot_output_path, "umap_refined_labels.png"),
             plot = p, width = 10, height = 8)
    },

    #' Plot a heatmap of custom gene list grouped by cluster
    #'
    #' This function plots a heatmap for a subset of genes, potentially grouped
    #' by cluster. It allows for customization of the color palette, scaling,
    #' and selection of top genes per cluster.
    #'
    #' @param genes                   Vector of genes to include in the heatmap.
    #' @param genes_by_cluster        Whether to group genes by cluster
    #'                                identity.
    #' @param n_top_genes_per_cluster Number of top genes per cluster to plot.
    #' @param color_palette           Custom color palette for clusters; uses
    #'                                hue_pal by default.
    #' @param scaled                  Whether the data is already scaled;
    #'                                applies
    #'                                row-wise z-score scaling if FALSE.
    #'
    #' @return                        Heatmap plot.
    plot_gene_heatmap = function(genes, genes_by_cluster = TRUE,
                                 n_top_genes_per_cluster = 5,
                                 color_palette = NULL, scaled = TRUE,
                                 assay = NULL, max_cells = 10000) {
      # assay was HARDCODED to "SCT". That is right in Step 4 and wrong in
      # Step 7, where find_top_genes returns VIPER proteins but the heatmap
      # looked their names up in SCT gene expression - so the figure showed
      # Pearson residuals of the TF genes, not NES. Two visible symptoms in job
      # 71416: mostly-blue with sparse red (TFs are lowly expressed, so most
      # cells sit at the residual floor), and an asymmetric ~-1..+4 colour scale
      # (Pearson residuals are bounded below and right-skewed, whereas z-scored
      # NES is roughly symmetric). It also silently dropped 17 proteins absent
      # from SCT - "Excluded genes: CCL4L2, PRDX2, ID3, ..." - which is the
      # signature of matching against the wrong assay, since every protein is
      # present in the VIPER assay by construction.
      #
      # Same class of bug as the hardcoded "integrated_snn_res." in Clusterer.
      if (is.null(assay)) assay <- DefaultAssay(self$seurat_obj)
      data_mat <-
        GetAssayData(self$seurat_obj, assay = assay, layer = "scale.data")
      if (is.null(data_mat) || nrow(data_mat) == 0) {
        stop("Assay '", assay, "' has no scale.data. Run ScaleData on it ",
             "before plotting, or pass a different `assay`.")
      }
      message(sprintf("  heatmap: assay '%s', %d x %d scale.data",
                      assay, nrow(data_mat), ncol(data_mat)))
      clusters <- self$seurat_obj$seurat_clusters

      if (length(unique(clusters)) == 0) {
        stop("No valid cluster data found.")
      }

      identities <- levels(factor(clusters))
      my_color_palette <-
        private$generate_color_palette(identities, color_palette)
      genes_in_data <- private$filter_genes(genes, data_mat)
      subset_mat <- private$subset_data(data_mat, genes_in_data, max_cells)
      cluster_data <- private$prepare_cluster_data(clusters, subset_mat)
      cluster_df <- cluster_data$cluster_df
      subset_mat <- cluster_data$subset_mat

      if (!scaled) {
        subset_mat <- private$apply_scaling(subset_mat)
      }

      mat_breaks <- private$generate_mat_breaks(subset_mat)
      annotations <-
        private$generate_annotations(cluster_df, my_color_palette,
                                     genes_by_cluster, n_top_genes_per_cluster)
      row_gaps <- private$adjust_row_gaps(annotations, cluster_df, subset_mat)

      pheatmap_args <-
        private$configure_pheatmap_args(subset_mat, cluster_df, mat_breaks,
                                        annotations, row_gaps, genes_by_cluster)
      heatmap_plot <- private$create_heatmap(pheatmap_args)
      private$save_heatmap(heatmap_plot, subset_mat)

      return(heatmap_plot)
    },

    #' Plot Cluster Frequencies by Provided Metadata
    #'
    #' @param plot_title       Title of the plot.
    #' @param group_by         Metadata label to group by. Any number of levels
    #'                         is supported.
    #' @param plot_type        Type of plot: "dot" for dot plot or "box" for
    #'                         box-whisker plot.
    #' @param binwidth         Width of bins in the dot plot.
    #' @param group_levels     Optional character vector fixing the order the
    #'                         groups appear in the legend and on the x axis.
    #'                         Defaults to the order in the data. Levels named
    #'                         here that are absent from the data are dropped
    #'                         with a warning rather than producing empty
    #'                         dodge slots.
    #'
    #' @return                 None. The function saves a plot to the specified
    #'                         directory.
    #'
    #' @note GENERALISED FROM TWO GROUPS TO N. This previously took the 1st and
    #'       2nd slices of the id x cluster x group table and combined them as
    #'       "early"/"late", so a third level was silently discarded - on this
    #'       project's neutrophil data, `type` is adjacentNormal/normal/tumor
    #'       and the 18 tumor Origins would have vanished from the figure with
    #'       no error. The `col_names` argument is also gone: it existed only
    #'       to smuggle the group label into a column name for melt() to
    #'       recover by splitting on "_", which broke on any id containing an
    #'       underscore and had to be kept in sync with the patient count by
    #'       hand. Group membership is now carried directly.
    plot_cluster_freq_by = function(plot_title, group_by,
                                    plot_type = "dot",
                                    binwidth = 0.01,
                                    group_levels = NULL) {
      required_columns <- c("id", "seurat_clusters", group_by)
      metadata_columns <-
        private$check_and_get_metadata_columns(self$seurat_obj,
                                               required_columns)

      id_column <- metadata_columns[["id"]]
      cluster_column <- metadata_columns[["seurat_clusters"]]
      group_by_column <- metadata_columns[[group_by]]

      cluster_freq_table <- table(id_column, cluster_column, group_by_column)
      plot_data <- private$cluster_freq_long(cluster_freq_table, group_levels)

      message(sprintf("Cluster frequency by '%s': %d group(s) [%s], %d id(s)",
                      group_by, nlevels(plot_data$type),
                      paste(levels(plot_data$type), collapse = ", "),
                      length(unique(plot_data$id))))

      output_file <-
        paste0("cluster_frequencies_by_", group_by, "_", plot_type, ".png")
      private$plot_cluster_frequencies(plot_data, plot_title,
                                       file.path(self$plot_output_path,
                                                 output_file),
                                       binwidth, plot_type)
    }
  ),

  #############################################################################
  #                           PRIVATE METHODS                                 #
  #############################################################################
  private = list(
    #' Filter Blueprint Labels Based on P-values and Frequency
    #'
    #' This function refines cell type labels based on p-values and frequency.
    #' Labels with p-values greater than 0.1 or occurring less than 50 times
    #' are set to NA.

    #' @return A Seurat object with refined labels.
    filter_blueprint_labels = function() {
      if (!all(c("blueprint_labels", "blueprint_pvals") %in%
                 colnames(self$seurat_obj@meta.data))) {
        stop(paste0("The Seurat object does not contain blueprint_labels",
                    " or blueprint_pvals. Please ensure these columns exist."))
      }

      refined_labels <- self$seurat_obj$blueprint_labels
      refined_labels[self$seurat_obj$blueprint_pvals > 0.1] <- NA
      refined_labels[refined_labels %in%
                       names(which(table(refined_labels) < 50))] <- NA
      self$seurat_obj$refined_labels <- refined_labels

      return(self$seurat_obj)
    },

    #' Generate a color palette
    #'
    #' Generates a color palette based on the provided identities. Defaults to
    #' hue_pal if no custom palette is provided.
    #'
    #' @param identities    Vector of unique identities for which colors are
    #'                      needed.
    #' @param color_palette Optional custom color palette.
    #'
    #' @return              A color palette vector.
    generate_color_palette = function(identities, color_palette = NULL) {
      if (is.null(color_palette)) {
        if (length(identities) > 0) {
          return(hue_pal()(length(identities)))
        } else {
          warning("No identities provided, returning empty color palette.")
          return(character(0))
        }
      } else {
        return(color_palette)
      }
    },

    #' Filter Genes Present in Data Matrix
    #'
    #' This function filters a given list of genes to include only those that
    #' are present in the specified data matrix. It provides informative
    #' messages about any genes that are not found in the data matrix.
    #'
    #' @param genes    A character vector of gene names to be included in the
    #'                 heatmap.
    #' @param data_mat A data matrix (e.g., expression matrix) with genes as
    #'                 row names.
    #'
    #' @return         A character vector of genes that are present in the
    #'                 data matrix.
    #' @throws         Error if none of the specified genes are present in the
    #'                 data matrix.
    filter_genes = function(genes, data_mat) {
      genes_in_data <- genes[genes %in% rownames(data_mat)]
      if (length(genes_in_data) == 0) {
        stop("None of the specified genes are present in the data.")
      }
      if (length(genes_in_data) < length(genes)) {
        warning("Some genes are not present in the data and will be excluded.")
        excluded_genes <- setdiff(genes, genes_in_data)
        message("Excluded genes: ", paste(excluded_genes, collapse = ", "))
      }
      return(genes_in_data)
    },

    #' Subset Data Matrix by Genes and Samples
    #'
    #' This function subsets the data matrix to include only the specified
    #' genes  and a random sample of columns (samples). It ensures that the
    #' dimensions of the subsetted data match the number of specified genes.
    #'
    #' @param data_mat      A data matrix (e.g., expression matrix) with genes
    #'                      as row names and samples as column names.
    #' @param genes_in_data A character vector of gene names that are present
    #'                      in the data matrix and should be included in the
    #'                      subset.
    #'
    #' @return              A subsetted data matrix containing only the
    #'                      specified genes and a random sample of up to 10,000
    #'                      columns (samples).
    #' @throws              Error if the number of rows in the subsetted data
    #'                      matrix does not match the number of specified
    #'                      genes.
    subset_data = function(data_mat, genes_in_data, max_cells = 10000,
                           seed = 1234) {
      # Seeded so the same object redraws the same figure. Unseeded, every run
      # picked a different 10,000 columns.
      set.seed(seed)
      # Columns are subsampled, as in scripts/archive/viper_global.R. At 49,142
      # (neutrophil) or 203,516 (macrophage) cells every column is well under a
      # pixel wide, so the figure shows aliasing rather than structure.
      i <-
        sample(seq_len(ncol(data_mat)), min(max_cells, ncol(data_mat)),
               replace = FALSE)
      subset_mat <- data_mat[genes_in_data, i]
      if (nrow(subset_mat) != length(genes_in_data)) {
        stop("Subset data dimensions do not match the number of genes.")
      }
      return(subset_mat)
    },

    #' Prepare Cluster Data and Subset Matrix
    #'
    #' This function prepares the cluster data frame and orders the subset
    #' matrix columns based on the cluster assignments. It ensures that the
    #' subset matrix columns are ordered according to their cluster identity.
    #'
    #' @param clusters   A factor or vector indicating the cluster assignments
    #'                   for each sample.
    #' @param subset_mat A subsetted data matrix with genes as row names and
    #'                   samples as column names.
    #'
    #' @return           A list containing:
    #'                    - cluster_df: A data frame with cluster assignments
    #'                                  for each sample, ordered by cluster
    #'                                  identity.
    #'                    - subset_mat: The subsetted data matrix with columns
    #'                                  ordered by cluster identity.
    prepare_cluster_data = function(clusters, subset_mat) {
      cluster_df <- data.frame(cluster = clusters[colnames(subset_mat)])
      rownames(cluster_df) <- colnames(subset_mat)
      ordered_indices <- order(cluster_df$cluster)
      subset_mat <- subset_mat[, ordered_indices]
      cluster_df <- cluster_df[ordered_indices, , drop = FALSE]
      return(list(cluster_df = cluster_df, subset_mat = subset_mat))
    },

    #' Apply Row-wise Z-score Scaling to Data Matrix
    #'
    #' This function applies z-score normalization across rows of a subsetted
    #' data matrix. It scales each gene (row) in the matrix by subtracting the
    #' mean and dividing by the standard deviation.
    #'
    #' @param subset_mat A subsetted data matrix with genes as row names and
    #'                   samples as column names.
    #'
    #' @return           A scaled data matrix where each gene (row) has been
    #'                   z-score normalized.
    apply_scaling = function(subset_mat) {
      return(t(apply(subset_mat, 1, private$calculate_z_score)))
    },

    #' Calculate row-wise z-score
    #'
    #' Applies z-score normalization across rows of a matrix.
    #'
    #' @param data_mat A numeric matrix.
    #' @return         Matrix with row-wise z-scores.
    calculate_z_score = function(data_mat) {
      return((data_mat - mean(data_mat)) / sd(data_mat))
    },

    #' Generate matrix breaks based on quantiles
    #'
    #' Determines breaks for the heatmap color scale based on quantiles,
    #' excluding extreme values.
    #'
    #' @param t Numeric matrix for which to determine breaks.
    #'
    #' @return  Vector of breaks for the heatmap color scale.
    generate_mat_breaks = function(t) {
      quantile_breaks <- function(xs, n = 10) {
        breaks <- quantile(xs, probs = seq(0, 1, length.out = n + 1))
        breaks[!duplicated(breaks)]
      }

      lower_breaks <- quantile_breaks(t[t < 0], n = 10)
      upper_breaks <- quantile_breaks(t[t > 0], n = 10)
      c(lower_breaks, 0, upper_breaks)[-c(1, length(lower_breaks),
                                          length(lower_breaks) + 2,
                                          length(lower_breaks) +
                                            length(upper_breaks) + 1)]
    },

    #' Generate annotations for heatmap
    #'
    #' Creates a list containing colors for cluster annotations and an optional
    #' data frame for row annotations if genes are grouped by cluster. The
    #' colors are matched to clusters, and if genes are grouped by cluster,
    #' each gene group is annotated with its corresponding cluster.
    #'
    #' @param df                      Data frame containing cluster information
    #'                                for columns in the heatmap.
    #' @param my_color_palette        Vector of colors used for cluster
    #'                                annotations.
    #' @param genes_by_cluster        Boolean indicating whether genes should
    #'                                be grouped by their cluster.
    #' @param n_top_genes_per_cluster Number of top genes per cluster to
    #'                                include if genes are grouped by cluster.
    #'
    #' @return                        A list containing 'anno_colors' for
    #'                                column annotations and 'anno_row' for row
    #'                                annotations (if applicable).
    generate_annotations = function(df, my_color_palette, genes_by_cluster,
                                    n_top_genes_per_cluster) {
      anno_colors <- list(cluster = my_color_palette)
      names(anno_colors$cluster) <- levels(df$cluster)

      if (genes_by_cluster) {
        anno_colors$group <- anno_colors$cluster
        anno_row <- data.frame(group = rep(levels(df$cluster),
                                           each = n_top_genes_per_cluster))
        return(list(anno_colors = anno_colors, anno_row = anno_row))
      } else {
        return(list(anno_colors = anno_colors, anno_row = NULL))
      }
    },

    #' Adjust Row Gaps for Heatmap Annotations
    #'
    #' This function calculates the row gaps for heatmap annotations based on
    #' the number of unique clusters and the actual number of top genes per
    #' cluster. It ensures that row gaps are correctly set for visual
    #' separation in the heatmap.
    #'
    #' @param annotations A list containing heatmap annotations, including row
    #'                    annotations.
    #' @param cluster_df  A data frame with cluster assignments for each
    #'                    sample, ordered by cluster identity.
    #' @param subset_mat  A subsetted data matrix with genes as row names and
    #'                    samples as column names.
    #'
    #' @return            A numeric vector of row gaps for heatmap
    #'                    visualization. Returns NULL if no row annotations are
    #'                    provided.
    adjust_row_gaps = function(annotations, cluster_df, subset_mat) {
      if (!is.null(annotations$anno_row)) {
        unique_clusters <- length(unique(cluster_df$cluster))
        n_top_genes_per_cluster_actual <-
          floor(nrow(subset_mat) / unique_clusters)
        row_gaps <- (2:unique_clusters - 1) * n_top_genes_per_cluster_actual
      } else {
        row_gaps <- NULL
      }
      return(row_gaps)
    },

    #' Configure pheatmap Arguments
    #'
    #' This function configures the arguments for the `pheatmap` function to
    #' create a heatmap. It sets various parameters including the data matrix,
    #' clustering options, annotations, and color settings.
    #'
    #' @param subset_mat       A subsetted data matrix with genes as row names
    #'                         and samples as column names.
    #' @param cluster_df       A data frame with cluster assignments for each
    #'                         sample, ordered by cluster identity.
    #' @param mat_breaks       A numeric vector of breakpoints for the heatmap
    #'                         color scale.
    #' @param annotations      A list containing heatmap annotations, including
    #'                         row and column annotations.
    #' @param gaps_row         A numeric vector of row gaps for heatmap
    #'                         visualization. Can be NULL.
    #' @param genes_by_cluster A boolean indicating whether to adjust the font
    #'                         size for row names based on clustering.
    #'
    #' @return                 A list of arguments configured for the
    #'                         `pheatmap` function.
    configure_pheatmap_args = function(subset_mat, cluster_df, mat_breaks,
                                       annotations, gaps_row,
                                       genes_by_cluster) {
      pheatmap_args <- list(subset_mat, cluster_rows = FALSE,
                            show_rownames = TRUE, cluster_cols = FALSE,
                            annotation_col = cluster_df, breaks = mat_breaks,
                            color = colorRampPalette(
                                                     c("blue",
                                                       "white",
                                                       "red"))
                            (length(mat_breaks)),
                            fontsize_row = ifelse(genes_by_cluster, 10, 8),
                            show_colnames = FALSE,
                            annotation_colors = annotations$anno_colors)

      if (!is.null(annotations$anno_row)) {
        pheatmap_args$annotation_row <- annotations$anno_row
        pheatmap_args$gaps_row <- gaps_row
      }

      return(pheatmap_args)
    },

    #' Create Heatmap Using pheatmap
    #'
    #' This function creates a heatmap using the `pheatmap` function with the
    #' specified arguments.
    #'
    #' @param pheatmap_args A list of arguments configured for the `pheatmap`
    #'                      function.
    #'
    #' @return              An object created by the `pheatmap` function, which
    #'                      contains the heatmap and associated metadata.
    create_heatmap = function(pheatmap_args) {
      return(do.call(pheatmap, pheatmap_args))
    },

    #' Save Heatmap to File
    #'
    #' This function saves the generated heatmap plot to a file with dynamically
    #' calculated dimensions.
    #'
    #' @param heatmap_plot An object created by the `pheatmap` function, which
    #'                     contains the heatmap and associated metadata.
    #' @param subset_mat   A subsetted data matrix with genes as row names and
    #'                     samples as column names.
    #'
    #' @return             None. The heatmap is saved to the specified file
    #'                     path.
    save_heatmap = function(heatmap_plot, subset_mat) {
      plot_width <- max(5, min(ncol(subset_mat) / 200, 25))
      plot_height <- max(8, min(nrow(subset_mat) / 5, 50))
      heatmap_plot_path <- file.path(self$plot_output_path, "gene_heatmap.png")
      ggsave(heatmap_plot_path, plot = heatmap_plot$gtable, width = plot_width,
             height = plot_height, limitsize = FALSE)
    },

    #' Filter and Subset Metadata Columns
    #'
    #' This function checks if the specified metadata columns are present in
    #' the Seurat object and returns the columns as a named list.
    #'
    #' @param seurat_obj   A Seurat object containing metadata columns.
    #' @param column_names A character vector of metadata column names.
    #'
    #' @return             A named list of metadata columns.
    check_and_get_metadata_columns = function(seurat_obj, column_names) {
      missing_columns <-
        column_names[!column_names %in% colnames(seurat_obj@meta.data)]
      if (length(missing_columns) > 0) {
        stop(paste("The following columns are missing in ",
                   "the Seurat object's metadata:",
                   paste(missing_columns, collapse = ", ")))
      }

      metadata_columns <-
        lapply(column_names, function(col) seurat_obj@meta.data[[col]])
      names(metadata_columns) <- column_names
      return(metadata_columns)
    },

    #' Calculate Cluster Frequencies
    #'
    #' @param cluster_freq_table Table of cluster frequencies.
    #' @param treatment_index Index of the treatment type.
    #'
    #' @return A data frame of normalized cluster frequencies.
    #' @note drop = FALSE is required. With a single id contributing to a group
    #'       the subset would collapse to a vector and apply() would then fail
    #'       or silently transpose. Groups contributed by one id are normal
    #'       here - several Origins are the only member of their Status.
    calculate_cluster_frequencies = function(cluster_freq_table,
                                             treatment_index) {
      freq_data <- as.data.frame.matrix(
        cluster_freq_table[, , treatment_index, drop = FALSE][, , 1])
      freq_data <- freq_data[rowSums(freq_data) > 0, , drop = FALSE]
      if (nrow(freq_data) == 0) {
        return(NULL)
      }
      # apply over rows returns clusters x ids (apply transposes).
      freq_data <- apply(freq_data, 1, function(x) {
        x / sum(x)
      })
      if (is.null(dim(freq_data))) {
        freq_data <- matrix(freq_data, ncol = 1,
                            dimnames = list(colnames(
                              cluster_freq_table[, , treatment_index,
                                                 drop = FALSE][, , 1]),
                              rownames(freq_data)))
      }
      return(freq_data)
    },

    #' Build a Long Data Frame of Per-Id Cluster Frequencies Across All Groups
    #'
    #' Replaces the old combine_cluster_frequencies + melt_cluster_frequencies
    #' pair, which only handled two groups and encoded group membership in
    #' column names.
    #'
    #' @param cluster_freq_table A 3-way table of id x cluster x group.
    #' @param group_levels       Optional ordering for the group factor.
    #'
    #' @return                   Data frame with columns cluster, type,
    #'                           frequency, id. One row per id per cluster per
    #'                           group the id contributes to.
    cluster_freq_long = function(cluster_freq_table, group_levels = NULL) {
      levels_present <- dimnames(cluster_freq_table)[[3]]
      if (is.null(levels_present)) {
        stop("The grouping variable produced no levels.")
      }

      pieces <- lapply(seq_along(levels_present), function(k) {
        freq <- private$calculate_cluster_frequencies(cluster_freq_table, k)
        if (is.null(freq)) {
          return(NULL)
        }
        # freq is clusters x ids
        data.frame(
          cluster = rep(rownames(freq), times = ncol(freq)),
          type = levels_present[k],
          frequency = as.vector(freq),
          id = rep(colnames(freq), each = nrow(freq)),
          stringsAsFactors = FALSE
        )
      })

      empty <- levels_present[vapply(pieces, is.null, logical(1))]
      if (length(empty) > 0) {
        warning("No cells for group(s): ", paste(empty, collapse = ", "))
      }
      out <- do.call(rbind, pieces)
      if (is.null(out) || nrow(out) == 0) {
        stop("No cluster frequencies could be computed for any group.")
      }

      keep <- intersect(if (is.null(group_levels)) levels_present
                        else group_levels,
                        unique(out$type))
      if (!is.null(group_levels)) {
        unknown <- setdiff(group_levels, unique(out$type))
        if (length(unknown) > 0) {
          warning("group_levels not present in the data, dropped: ",
                  paste(unknown, collapse = ", "))
        }
      }
      out <- out[out$type %in% keep, , drop = FALSE]
      out$type <- factor(out$type, levels = keep)
      out$cluster <- factor(out$cluster,
                            levels = dimnames(cluster_freq_table)[[2]])
      return(out)
    },

    #' Combine Cluster Frequencies
    #'
    #' @param early_data Data frame of early treatment cluster frequencies.
    #' @param late_data Data frame of late treatment cluster frequencies.
    #' @param col_names Vector of column names for the combined data frame.
    #'
    #' @return A combined data frame of early and late treatment cluster
    #'         frequencies.
    combine_cluster_frequencies = function(early_data, late_data, col_names) {
      combined <- cbind(early_data, late_data)
      colnames(combined) <- col_names
      return(combined)
    },

    #' Melt Cluster Frequencies
    #'
    #' @param combined_data Combined data frame of cluster frequencies.
    #'
    #' @return A melted data frame suitable for plotting.
    melt_cluster_frequencies = function(combined_data) {
      melted_data <- melt(combined_data)
      colnames(melted_data) <- c("cluster", "type", "frequency")
      melted_data$type <- unlist(
        lapply(strsplit(as.character(melted_data$type), "_"), function(x) {
          x[1]
        })
      )
      melted_data$type <- factor(melted_data$type, levels = c("Early", "Late"))
      melted_data$cluster <- as.factor(melted_data$cluster)
      return(melted_data)
    },

    #' Plot Cluster Frequencies
    #'
    #' @param plot_data   Data frame of cluster frequencies to be plotted.
    #' @param plot_title  Title of the plot.
    #' @param output_path File path to save the plot.
    #' @param binwidth    Width of bins in the dot plot.
    #' @param plot_type   Type of plot: "dot" for dot plot or "box" for
    #'                    box-whisker plot.
    #'
    #' @return            None. The function saves a plot to the specified file
    #'                    path.
    plot_cluster_frequencies =
      function(plot_data, plot_title, output_path, binwidth, plot_type) {
        theme_update(plot.title = element_text(hjust = 0.5))

        if (plot_type == "dot") {
          private$plot_dot_plot(plot_data, plot_title, output_path, binwidth)
        } else if (plot_type == "box") {
          private$plot_box_plot(plot_data, plot_title, output_path)
        } else {
          stop("Invalid plot type specified. Use 'dot' or 'box'.")
        }
      },

    #' Plot Dot Plot
    #'
    #' @param plot_data   Data frame of cluster frequencies to be plotted.
    #' @param plot_title  Title of the plot.
    #' @param output_path File path to save the plot.
    #' @param binwidth    Width of bins in the dot plot.
    #'
    #' @return            None. The function saves a dot plot to the specified
    #'                    file path.
    plot_dot_plot = function(plot_data, plot_title, output_path, binwidth) {
      p <- ggplot(plot_data, aes(x = cluster, y = frequency, fill = type)) +
        geom_dotplot(binaxis = "y", stackdir = "center",
                     position = position_dodge(), binwidth = binwidth) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")) +
        ggtitle(plot_title) +
        theme(plot.title = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 14, face = "bold"),
              axis.text = element_text(size = 8),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 12),
              strip.text.x = element_text(size = 12, face = "bold"))

      ggsave(output_path, plot = p, width = 10, height = 8)
    },
    #' Plot Box Plot
    #'
    #' @param plot_data   Data frame of cluster frequencies to be plotted.
    #' @param plot_title  Title of the plot.
    #' @param output_path File path to save the plot.
    #'
    #' @return            None. The function saves a box plot to the specified
    #'                    file path.
    plot_box_plot = function(plot_data, plot_title, output_path) {
      p <- ggplot(plot_data, aes(x = cluster, y = frequency, fill = type)) +
        geom_boxplot(position = position_dodge()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")) +
        ggtitle(plot_title) +
        theme(plot.title = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 14, face = "bold"),
              axis.text = element_text(size = 8),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 12),
              strip.text.x = element_text(size = 12, face = "bold"))

      ggsave(output_path, plot = p, width = 10, height = 8)
    }
  )
)