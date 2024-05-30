
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ========================================================
# Step 7: Re-clustering based on VIPER results
# ========================================================

# Verify that viper_results is a list and correctly formatted
if (!is.list(viper_results)) {
  stop("Error: VIPER results must be a list format.")
}

# Combine the list elements of viper_results into a matrix
viper_data <- do.call(cbind, viper_results)
colnames(viper_data) <- colnames(integrated_seurat)

# Output dimensions of viper_data to verify alignment with expected structure
cat("Dimensions of combined VIPER data matrix:")
print(dim(viper_data))

# The # of columns in viper_data matches the # of cells in integrated_seurat?
if (ncol(viper_data) != ncol(integrated_seurat)) {
  stop("Error: Mismatch in # of cells between VIPER and integrated data.",
       " Expected", ncol(integrated_seurat), "cells, but got",
       ncol(viper_data), "in VIPER results.")
}

# Validate that the column names in viper_data and integrated_seurat match
if (!all(colnames(viper_data) == colnames(integrated_seurat))) {
  stop(paste0("Error: Column names of the VIPER data do not",
              " match those of the Seurat object."))
}

# Add VIPER scores as a new assay in the Seurat object
integrated_seurat[["VIPER_scores"]] <- CreateAssayObject(viper_data)
DefaultAssay(integrated_seurat) <- "VIPER_scores"

feature_variances <- apply(viper_data, 1, var)

# Filter out features with zero variance
non_constant_features <- viper_data[feature_variances != 0, ]

# Debug: Print the number of features with non-zero variance
cat("Number of features with non-zero variance:",
    nrow(non_constant_features), "\n")

# Ensure there are enough features left to continue the analysis
if (nrow(non_constant_features) < 2) {
  stop("Insufficient variable features for PCA.")
}

# Update the assay object with non-constant features
integrated_seurat[["VIPER_scores"]] <- CreateAssayObject(non_constant_features)
DefaultAssay(integrated_seurat) <- "VIPER_scores"

# Perform hierarchical clustering on transposed data
dist_matrix <- dist(t(non_constant_features))
hc <- hclust(dist_matrix)
clusters <- cutree(hc, k = 5)

# Add cluster assignments to Seurat object
integrated_seurat <- AddMetaData(integrated_seurat, metadata = clusters,
                                 col.name = "seurat_clusters")

# Visualize the re-clustering results using UMAP
p <- DimPlot(integrated_seurat, reduction = "umap",
             group.by = "seurat_clusters") +
  ggtitle("Re-clustering Based on VIPER Scores")

# Define the path and filename for saving the UMAP plot of re-clustered data
reclustered_umap_plot_path <-
  file.path(plot_output_path, "reclustered_umap_results.png")
ggsave(reclustered_umap_plot_path, plot = p, width = 10, height = 8)

# ========================================================
# Step 8: Plotting the frequency of each cluster by patient
# ========================================================

plot_cluster_frequency <- function(data, cluster_label, patient_label) {
  df <- data@meta.data %>%
    dplyr::select({{cluster_label}}, {{patient_label}}) %>%
    dplyr::group_by(.data[[cluster_label]], .data[[patient_label]]) %>%
    dplyr::summarise(count = n(), .groups = "drop")

  ggplot(df,
         aes(x = as.factor(.data[[cluster_label]]), y = count,
             fill = as.factor(.data[[patient_label]]))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Cluster", y = "Count", fill = "Patient") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      panel.border = element_blank(), # Remove the panel border
      plot.background = element_rect(fill = "white", colour = "white")
    ) +
    ggtitle("Frequency of Each Cluster by Patient")
}

# Function call to plot cluster frequency
p1 <- plot_cluster_frequency(integrated_seurat,
                             cluster_label = "seurat_clusters",
                             patient_label = "patient")

# Define the path and filename for saving the cluster frequency plot
cluster_frequency_plot_path <-
  file.path(plot_output_path, "cluster_frequency_by_patient.png")
ggsave(cluster_frequency_plot_path, plot = p1, width = 10, height = 8)

# Save the final integrated and annotated Seurat object
saveRDS(integrated_seurat,
        file = file.path(base_output_path,
                         "final_integrated_seurat_object.rds"))

# Identify top genes per cluster if not predefined
top_genes_per_cluster <- FindAllMarkers(integrated_seurat, only.pos = TRUE,
                                        min.pct = 0.25,
                                        logfc.threshold = 0.25)

# Check the availability of genes per cluster
gene_counts <- top_genes_per_cluster %>%
  group_by(cluster) %>%
  summarise(n_genes = n())

# Filter out clusters with fewer than required genes
valid_clusters <- gene_counts %>% filter(n_genes >= 5)

# Filter top_genes to include only those in valid clusters
top_genes <- top_genes_per_cluster %>%
  filter(cluster %in% valid_clusters$cluster) %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Ensure all selected clusters have enough genes
if (nrow(top_genes) < length(unique(valid_clusters$cluster)) * 5) {
  cat("Not all clusters have enough top genes for the heatmap.\n")
} else {
  # Check for duplicates in the gene list
  if (length(unique(top_genes$gene)) != length(top_genes$gene)) {
    cat("Duplicate gene names found in the top genes list.\n")
  }

  # Attempt the heatmap plot if the gene list is correct
  gene_data <- GetAssayData(integrated_seurat, layer = "data")
  if (any(!top_genes$gene %in% rownames(gene_data))) {
    cat("Some genes in top_genes not found in the data matrix.\n")
  } else {
    p2 <- gene_heatmap_plot(gene_data,
                            integrated_seurat$seurat_clusters,
                            genes = top_genes$gene,
                            n_top_genes_per_cluster = 5,
                            scaled = TRUE) +
      ggtitle("Heatmap of Top 5 Genes Per Cluster")

    # Define the path and filename for saving the gene heatmap plot
    gene_heatmap_plot_path <-
      file.path(plot_output_path, "gene_heatmap_per_cluster.png")
    ggsave(gene_heatmap_plot_path, plot = p2, width = 10, height = 10)
  }
}