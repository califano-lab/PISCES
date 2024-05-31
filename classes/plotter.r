library(R6)
library(Seurat)
library(ggplot2)

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
      errbar(self$resolutions, mean_scores, mean_scores + sd_scores,
             mean_scores - sd_scores, ylab = "Mean Silhouette Score",
             xlab = "Resolution Parameter")
      lines(self$resolutions, mean_scores)

      best_resolution <-
        tail(self$resolutions[which(mean_scores == max(mean_scores))], n = 1)
      legend("topright",
             paste("Best Resolution", best_resolution, sep = " = "))

      plot_path <- file.path(self$plot_output_path, "silhouette_scores.png")
      dev.copy(png, filename = plot_path)
      dev.off()
    },

    #' Plot UMAP Clusters
    #'
    #' This function plots UMAP clusters for a Seurat object, assigning
    #' provided cluster labels and saving the plot to the specified output
    #' path.
    #'
    #' @return The Seurat object with updated cluster labels.
    plot_umap_clusters = function() {
      unique_clusters <- unique(self$seurat_obj$seurat_clusters)
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

      # Map numeric cluster IDs to the provided labels
      self$seurat_obj$seurat_clusters <-
        mapvalues(self$seurat_obj$seurat_clusters,
                  from = seq_along(unique_clusters) - 1,
                  to = self$cluster_labels)

      # Create and save the UMAP plot
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
    #' @todo                          Ask doctor of the gene exclusion and
    #'                                and refactor further.
    plot_gene_heatmap = function(genes, genes_by_cluster = TRUE,
                                 n_top_genes_per_cluster = 5,
                                 color_palette = NULL, scaled = FALSE) {
      dat <- GetAssayData(self$seurat_obj, assay = "SCT", layer = "scale.data")
      clust <- self$seurat_obj$seurat_clusters

      if (length(unique(clust)) == 0) {
        stop("No valid cluster data found.")
      }
      identities <- levels(factor(clust))

      # Prepare color palette
      my_color_palette <-
        private$generate_color_palette(identities, color_palette)

      # Filter genes to include only those present in the data
      genes_in_data <- genes[genes %in% rownames(dat)]
      if (length(genes_in_data) == 0) {
        stop("None of the specified genes are present in the data.")
      }
      if (length(genes_in_data) < length(genes)) {
        warning("Some genes are not present in the data and will be excluded.")
        excluded_genes <- setdiff(genes, genes_in_data)
        message("Excluded genes: ", paste(excluded_genes, collapse = ", "))
      }

      # Subset data for heatmap
      i <- sample(seq_len(ncol(dat)), min(10000, ncol(dat)), replace = FALSE)
      x <- dat[genes_in_data, i]

      # Validate dimensions after subsetting
      if (nrow(x) != length(genes_in_data)) {
        stop("Subset data dimensions do not match the number of genes.")
      }

      # Prepare cluster data frame
      df <- data.frame(cluster = clust[i])
      rownames(df) <- colnames(x)
      o <- order(df$cluster)
      x <- x[, o]
      df <- df[o, , drop = FALSE]

      # Apply scaling if needed
      if (!scaled) {
        x <- t(apply(x, 1, private$calculate_z_score))
      }

      # Generate breaks and annotations
      mat_breaks <- private$generate_mat_breaks(x)
      annotations <-
        private$generate_annotations(df, my_color_palette, genes_by_cluster,
                                     n_top_genes_per_cluster)

      # Adjust gaps_row to match the actual number of genes
      if (!is.null(annotations$anno_row)) {
        unique_clusters <- length(unique(df$cluster))
        n_top_genes_per_cluster_actual <- floor(nrow(x) / unique_clusters)
        gaps_row <- (2:unique_clusters - 1) * n_top_genes_per_cluster_actual
      } else {
        gaps_row <- NULL
      }

      # Configure pheatmap arguments
      pheatmap_args <- list(x, cluster_rows = FALSE, show_rownames = TRUE,
                            cluster_cols = FALSE, annotation_col = df,
                            breaks = mat_breaks,
                            color = colorRampPalette(c("blue",
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

      # Create heatmap
      heatmap_plot <- do.call(pheatmap, pheatmap_args)

      # Save heatmap to file
      heatmap_plot_path <- file.path(self$plot_output_path, "gene_heatmap.png")
      ggsave(heatmap_plot_path, plot = heatmap_plot$gtable, width = 10,
             height = 8)

      return(heatmap_plot)
    },

    #' Plot Cluster Frequencies by Treatment
    #'
    #' @param col_names        Vector of column names for the combined data
    #'                         frame.
    #' @param plot_title       Title of the plot.
    #' @param binwidth         Width of bins in the dot plot.
    #'
    #' @return                 None. The function saves a plot to the specified
    #'                         directory.
    #' @todo                   Abstract so that it can be used for other...
    plot_cluster_freq_by_treatment = function(col_names, plot_title,
                                              binwidth = 0.01) {

      cluster_freq_table <-
        table(self$seurat_obj$id, self$seurat_obj$seurat_clusters,
              self$seurat_obj$type)

      early_data <-
        private$calculate_cluster_frequencies(cluster_freq_table, 1)
      late_data <-
        private$calculate_cluster_frequencies(cluster_freq_table, 2)

      combined_data <-
        private$combine_cluster_frequencies(early_data, late_data, col_names)

      plot_data <- private$melt_cluster_frequencies(combined_data)

      private$plot_cluster_frequencies(plot_data, plot_title,
                                       file.path(
                                         self$plot_output_path,
                                         "cluster_frequencies_by_treatment.png"
                                       ),
                                       binwidth)
    }
  ),

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

    #' Calculate row-wise z-score
    #'
    #' Applies z-score normalization across rows of a matrix.
    #'
    #' @param x A numeric matrix.
    #' @return  Matrix with row-wise z-scores.
    calculate_z_score = function(x) {
      return((x - mean(x)) / sd(x))
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

    #' Calculate Cluster Frequencies
    #'
    #' @param cluster_freq_table Table of cluster frequencies.
    #' @param treatment_index Index of the treatment type.
    #'
    #' @return A data frame of normalized cluster frequencies.
    calculate_cluster_frequencies = function(cluster_freq_table,
                                             treatment_index) {
      freq_data <- as.data.frame.matrix(cluster_freq_table[, , treatment_index])
      freq_data <- freq_data[rowSums(freq_data) > 0, ]
      freq_data <- apply(freq_data, 1, function(x) {
        x / sum(x)
      })
      return(freq_data)
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
    #' @param plot_data Data frame of cluster frequencies to be plotted.
    #' @param plot_title Title of the plot.
    #' @param output_path File path to save the plot.
    #' @param binwidth Width of bins in the dot plot.
    #'
    #' @return None. The function saves a plot to the specified file path.
    plot_cluster_frequencies =
      function(plot_data, plot_title, output_path, binwidth) {
        theme_update(plot.title = element_text(hjust = 0.5))

        p <- ggplot(plot_data, aes(x = cluster, y = frequency, fill = type)) + # nolint
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
      }
  )
)