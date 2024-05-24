library(dplyr)
library(cluster)
library(umap)
library(pheatmap)
library(Hmisc)
library(Seurat)
library(plyr)
library(ggplot2)
library(scales)

#' Find the Best Clustering Resolution Using Silhouette Scores
#'
#' This function identifies the best clustering resolution for a Seurat object
#' by calculating silhouette scores for various resolutions and selecting the
#' one with the highest mean silhouette score.
#'
#' @param seurat_obj A Seurat object containing the clustering results and PCA
#'                   embeddings.
#' @param resolutions A vector of clustering resolutions to evaluate.
#' @param pca_dims A vector of PCA dimensions to use for silhouette score
#'                 calculation.
#'
#' @return A list containing:
#'  - best_resolution: The resolution with the highest mean silhouette score.
#'  - mean_scores: A named vector of mean silhouette scores for each
#'                 resolution.
#'  - sd_scores: A named vector of standard deviations of silhouette scores for
#'               each resolution.
find_best_resolution <- function(seurat_obj, resolutions, pca_dims) {
  # Extract clustering results for all resolutions
  clustering_results <-
    seurat_obj@meta.data[, grepl("integrated_snn_res.",
                                 colnames(seurat_obj@meta.data))]

  # Extract PCA embeddings as a matrix
  pca_matrix <-
    as.data.frame(t(seurat_obj@reductions$pca@cell.embeddings[pca_dims, ]))

  # Calculate silhouette scores for each resolution
  silhouette_scores <-
    compute_silhouette_scores(pca_matrix, clustering_results)
  mean_scores <- silhouette_scores$means
  sd_scores <- silhouette_scores$sd

  # Identify the best resolution based on the highest mean silhouette score
  best_resolution <-
    tail(resolutions[which(mean_scores == max(mean_scores))], n = 1)

  return(list(best_resolution = best_resolution, mean_scores = mean_scores,
              sd_scores = sd_scores))
}

#' Compute silhouette scores for multiple Louvain clusterings
#'
#' This function evaluates multiple Louvain clusterings of a single-cell matrix
#' by computing the silhouette scores for a subsample of the data. It iterates
#' over 100 alternative resolution values, each time sub-sampling 1000 cells
#' (or fewer if fewer are available) and calculates the mean and standard
#' deviation of the silhouette scores for each clustering resolution.
#'
#' @param mat   Matrix with rows as principal component vectors and columns as
#'              samples.
#' @param clust Matrix with rows as samples and each column as a clustering
#'              vector for a given resolution.
#'
#' @return      List containing the means and standard deviations of silhouette
#'              scores for each clustering resolution.
compute_silhouette_scores <- function(mat, clust) {
  num_resolutions <- ncol(clust)
  num_subsamples <- 100

  # Initialize output dataframe
  silhouette_scores <- initialize_silhouette_scores(num_subsamples,
                                                    num_resolutions)

  for (resolution_index in 1:num_resolutions) {
    for (subsample_index in 1:num_subsamples) {
      sampled_indices <- sample_cells(mat, 1000)
      distance_matrix <- compute_distance_matrix(mat[, sampled_indices])
      silhouette_scores[subsample_index, resolution_index] <-
        compute_silhouette_width(clust[sampled_indices, resolution_index],
                                 distance_matrix)
    }
  }

  list(means = colMeans(silhouette_scores, na.rm = TRUE),
       sd = apply(silhouette_scores, 2, sd, na.rm = TRUE))
}

#' Initialize matrix for storing silhouette scores
#'
#' @param num_subsamples  Number of subsamples to compute.
#' @param num_resolutions Number of resolution values/clustering vectors.
#'
#' @return                Initialized matrix for storing silhouette scores.
initialize_silhouette_scores <- function(num_subsamples, num_resolutions) {
  matrix(rep(NA, num_subsamples * num_resolutions), nrow = num_subsamples)
}

#' Randomly sample cells from the matrix
#'
#' @param mat       Data matrix.
#' @param num_cells Number of cells to sample.
#'
#' @return          Indices of sampled cells.
sample_cells <- function(mat, num_cells) {
  sample(seq_len(ncol(mat)), min(num_cells, ncol(mat)))
}

#' Compute Distance Matrix
#'
#' Calculates a distance matrix using Pearson correlation.
#'
#' @param dat_mat A matrix of gene expression data (genes x samples).
#'
#' @return        A distance matrix.
compute_distance_matrix <- function(dat_mat) {
  if (!is.matrix(dat_mat)) {
    dat_mat <- as.matrix(dat_mat)
  }

  dist_mat <- as.dist(1 - cor(dat_mat, method = "pearson"))
  return(dist_mat)
}

#' Compute the mean silhouette width for a clustering
#'
#' @param clustering      Clustering vector for the subsampled cells.
#' @param distance_matrix Distance matrix for the subsampled cells.
#'
#' @return                Mean silhouette width for the clustering.
compute_silhouette_width <- function(clustering, distance_matrix) {
  if (length(unique(clustering)) <= 1) return(0)

  silhouette_scores <- silhouette(as.numeric(clustering), distance_matrix)
  mean(silhouette_scores[, "sil_width"])
}

#' Plot Silhouette Scores
#'
#' This function plots the mean silhouette scores with error bars for different
#' clustering resolutions and highlights the best resolution based on the
#' highest mean silhouette score.
#'
#' @param resolutions A vector of clustering resolutions.
#' @param mean_scores A vector of mean silhouette scores for each resolution.
#' @param sd_scores A vector of standard deviations of silhouette scores for
#'                  each resolution.
#' @param plot_output_path The directory path where the plot will be saved.
#'
#' @return None. The function saves the plot as a PNG file in the specified
#'         directory.
plot_silhouette_scores <-
  function(resolutions, mean_scores, sd_scores, plot_output_path) {
    errbar(resolutions, mean_scores, mean_scores + sd_scores,
           mean_scores - sd_scores, ylab = "Mean Silhouette Score",
           xlab = "Resolution Parameter")
    lines(resolutions, mean_scores)

    best_resolution <-
      tail(resolutions[which(mean_scores == max(mean_scores))], n = 1)
    legend("topright", paste("Best Resolution", best_resolution, sep = " = "))

    plot_path <- file.path(plot_output_path, "silhouette_scores.png")
    dev.copy(png, filename = plot_path)
    dev.off()
  }

#' Set Clusters Based on the Best Resolution
#'
#' This function assigns clusters to a Seurat object based on the best
#' clustering resolution and sets the active identity class to these clusters.
#'
#' @param seurat_obj A Seurat object containing the clustering results.
#' @param best_resolution The best clustering resolution identified.
#'
#' @return The Seurat object with clusters set based on the best resolution.
set_best_clusters <- function(seurat_obj, best_resolution) {
  # Assign clusters based on the best resolution
  cluster_column_name <- paste("integrated_snn_res.", best_resolution, sep = "")
  seurat_obj$seurat_clusters <- seurat_obj@meta.data[, cluster_column_name]

  # Set the active identity class to the new clusters
  Idents(seurat_obj) <- "seurat_clusters"

  return(seurat_obj)
}

#' Plot UMAP Clusters
#'
#' This function plots UMAP clusters for a Seurat object, assigning provided
#' cluster labels and saving the plot to the specified output path.
#'
#' @param seurat_obj A Seurat object containing UMAP and clustering results.
#' @param cluster_labels A character vector of cluster labels to assign to the
#'                       clusters.
#' @param plot_output_path The path where the UMAP plot will be saved.
#'
#' @return The Seurat object with updated cluster labels.
plot_umap_clusters <- function(seurat_obj, cluster_labels, plot_output_path) {
  unique_clusters <- unique(seurat_obj$seurat_clusters)
  num_clusters <- length(unique_clusters)

  # Ensure the number of provided labels matches the number of clusters
  if (num_clusters > length(cluster_labels)) {
    warning(paste0("More clusters found in data than provided labels. ",
                   "Some clusters will be unlabeled."))
    cluster_labels <-
      c(cluster_labels,
        paste0("Cluster", (length(cluster_labels) + 1):num_clusters))
  } else if (num_clusters < length(cluster_labels)) {
    cluster_labels <- cluster_labels[1:num_clusters]
  }

  # Map numeric cluster IDs to the provided labels
  seurat_obj$seurat_clusters <- mapvalues(seurat_obj$seurat_clusters,
                                          from = seq_along(unique_clusters) - 1,
                                          to = cluster_labels)

  # Create and save the UMAP plot
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters",
               label = TRUE, label.size = 7, repel = TRUE) + NoLegend()
  umap_cluster_plot_path <-
    file.path(plot_output_path, "umap_clustering_results.png")
  ggsave(umap_cluster_plot_path, plot = p, width = 10, height = 8)

  return(seurat_obj)
}

#' Find top genes for each cluster
#'
#' This function identifies the top marker genes for each cluster in the
#' integrated Seurat object.
#'
#' @param seurat_obj A Seurat object with clustering information.
#' @param assay_name The name of the assay to use for finding markers
#' @param n_top_genes The number of top genes to find for each cluster
#' @param logfc_threshold The log fold change threshold for marker genes
#'
#' @return A data frame containing the top marker genes for each cluster.
find_top_genes <- function(seurat_obj, assay_name = "SCT", n_top_genes = 10,
                           logfc_threshold = 0.25) {
  # Prepare the SCT assay for differential expression analysis
  seurat_obj <- PrepSCTFindMarkers(seurat_obj)

  # Find all markers
  all_markers <- FindAllMarkers(
    seurat_obj,
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

#' Plot a heatmap of custom gene list grouped by cluster
#'
#' This function plots a heatmap for a subset of genes, potentially grouped by
#' cluster. It allows for customization of the color palette, scaling, and
#' selection of top genes per cluster.
#'
#' @param dat                     Matrix with genes as rows and samples as
#'                                columns.
#' @param clust                   Vector of cluster labels.
#' @param genes                   Vector of genes to include in the heatmap.
#' @param genes_by_cluster        Whether to group genes by cluster identity.
#' @param n_top_genes_per_cluster Number of top genes per cluster to plot.
#' @param color_palette           Custom color palette for clusters; uses
#'                                hue_pal by default.
#' @param scaled                  Whether the data is already scaled; applies
#'                                row-wise z-score scaling if FALSE.
#' @param plot_output_path        Path to save the heatmap plot.
#'
#' @return                        Heatmap plot.
#' @todo                          Ask doctor of the gene exclusion and
#'                                and refactor furhter
plot_gene_heatmap <- function(dat, clust, genes, genes_by_cluster = TRUE,
                              n_top_genes_per_cluster = 5, color_palette = NULL,
                              scaled = FALSE, plot_output_path) {

  if (length(unique(clust)) == 0) {
    stop("No valid cluster data found.")
  }
  identities <- levels(factor(clust))

  # Prepare color palette
  my_color_palette <- generate_color_palette(identities, color_palette)

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
    x <- t(apply(x, 1, calculate_z_score))
  }

  # Generate breaks and annotations
  mat_breaks <- generate_mat_breaks(x)
  annotations <- generate_annotations(df, my_color_palette, genes_by_cluster,
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
  heatmap_plot_path <- file.path(plot_output_path, "gene_heatmap.png")
  ggsave(heatmap_plot_path, plot = heatmap_plot$gtable, width = 10, height = 8)

  return(heatmap_plot)
}

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
generate_color_palette <- function(identities, color_palette = NULL) {
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
}

#' Calculate row-wise z-score
#'
#' Applies z-score normalization across rows of a matrix.
#'
#' @param x A numeric matrix.
#' @return  Matrix with row-wise z-scores.
calculate_z_score <- function(x) {
  return((x - mean(x)) / sd(x))
}

#' Generate matrix breaks based on quantiles
#'
#' Determines breaks for the heatmap color scale based on quantiles, excluding
#' extreme values.
#'
#' @param t Numeric matrix for which to determine breaks.
#'
#' @return  Vector of breaks for the heatmap color scale.
generate_mat_breaks <- function(t) {
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
}

#' Generate annotations for heatmap
#'
#' Creates a list containing colors for cluster annotations and an optional
#' data frame for row annotations if genes are grouped by cluster. The colors
#' are matched to clusters, and if genes are grouped by cluster, each gene
#' group is annotated with its corresponding cluster.
#'
#' @param df                      Data frame containing cluster information for
#'                                columns in the heatmap.
#' @param my_color_palette        Vector of colors used for cluster
#'                                annotations.
#' @param genes_by_cluster        Boolean indicating whether genes should be
#'                                grouped by their cluster.
#' @param n_top_genes_per_cluster Number of top genes per cluster to include if
#'                                genes are grouped by cluster.
#'
#' @return                        A list containing 'anno_colors' for column
#'                                annotations and 'anno_row' for row
#'                                annotations (if applicable).
generate_annotations <- function(df, my_color_palette, genes_by_cluster,
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
}

#' Filter Blueprint Labels Based on P-values and Frequency
#'
#' This function refines cell type labels based on p-values and frequency.
#' Labels with p-values greater than 0.1 or occurring less than 50 times are
#' set to NA.
#'
#' @param seurat_obj A Seurat object with blueprint labels and p-values in the
#'                   metadata.
#' @return A Seurat object with refined labels.
#' @todo Ask doctor about strictness of p-value and frequency thresholds.
filter_blueprint_labels <- function(seurat_obj) {
  if (!all(c("blueprint_labels", "blueprint_pvals") %in%
             colnames(seurat_obj@meta.data))) {
    stop(paste0("The Seurat object does not contain blueprint_labels",
                " or blueprint_pvals. Please ensure these columns exist."))
  }

  refined_labels <- seurat_obj$blueprint_labels
  refined_labels[seurat_obj$blueprint_pvals > 0.1] <- NA
  refined_labels[refined_labels %in%
                   names(which(table(refined_labels) < 50))] <- NA
  seurat_obj$refined_labels <- refined_labels

  return(seurat_obj)
}

#' Plot UMAP with Refined Labels
#'
#' This function plots a UMAP visualization of the Seurat object with refined
#' labels.
#'
#' @param seurat_obj A Seurat object with refined labels in the metadata.
#' @param plot_output_path The directory path where the plot should be saved.
#' @return The UMAP plot is saved to the specified directory.
plot_umap_with_labels <- function(seurat_obj, plot_output_path) {
  p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE,
               label.size = 5, group.by = "refined_labels") + NoLegend()
  ggsave(file.path(plot_output_path, "umap_refined_labels.png"), plot = p,
         width = 10, height = 8)
}