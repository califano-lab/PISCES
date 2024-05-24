source("standard_workflow_code.r")

library(Atools)
library(MASS)
library(reticulate)
library(SingleR)
library(Seurat)
library(cowplot)
library(dplyr)
library(cluster)
library(umap)
library(reshape)
library(pheatmap)
library(viper)
library(clustree)
library(factoextra)
library(Hmisc)
library(ggplot2)
library(scales)
library(ggrepel)
library(plyr)

################## DEFINE YOUR LOCAL PATHS HERE ##################

base_path <- "/Users/apple/Desktop/240307_JOEL_DAVID_6_HUMAN_10X"
base_output_path <- "/Users/apple/Desktop/output"
plot_output_path <- file.path(base_output_path, "plots")
aracne_binary_path <- paste0("/Users/apple/Documents/Research/aleks-lab/",
                             "repos/ARACNe3/build/src/app/",
                             "ARACNe3_app_release")

# Define paths to regulator files
regulator_dir <- file.path(base_output_path, "human_hugo")
regulator_files <- list(
  cotfs = file.path(regulator_dir, "cotfs-hugo.txt"),
  surface = file.path(regulator_dir, "surface-hugo.txt"),
  sig = file.path(regulator_dir, "sig-hugo.txt"),
  tfs = file.path(regulator_dir, "tfs-hugo.txt")
)

# Create the directory if it does not exist
if (!dir.exists(plot_output_path)) {
  dir.create(plot_output_path, recursive = TRUE)
}
#################################################################

################## DEFINE OTHER PREFERENCES #####################

my_verbose <- FALSE

#################################################################

# ========================================================
# Step 1: Load and preprocess data
# ========================================================

# Define patient information
patients <- list(
  list(id = "JD001", type = "Early"),
  list(id = "JD002", type = "Early"),
  list(id = "JD003", type = "Late"),
  list(id = "JD004", type = "Late"),
  list(id = "JD005", type = "Early"),
  list(id = "JD006", type = "Late")
)

# Define constants for the data path construction
analysis_prefix <- "JOEL_DAVID_6_HUMAN_10X-"
count_default_suffix <- "-cellranger-count-default"
output_folder_suffix <- "_cellranger_count_outs"
feature_matrix_dir <- "filtered_feature_bc_matrix"

# Load data for each patient into separate Seurat objects
patient_seurat_list <- lapply(patients, function(patient) {
  tryCatch({
    message("Loading patient: ", patient$id)
    seurat_obj <- load_into_seurat(patient, base_path, analysis_prefix,
                                   count_default_suffix, output_folder_suffix,
                                   feature_matrix_dir)
    message("Completed loading for patient: ", patient$id)
    return(seurat_obj)
  }, error = function(e) {
    stop("Error loading patient: ", patient$id, ": ", e$message)
  })
})

message("Finished loading patient data into Seurat objects.")

blueprint_encode <- BlueprintEncodeData()

patient_seurat_list <- lapply(patient_seurat_list, function(p) {
  patient_id <- unique(p$patient)
  tryCatch({
    message("Preprocessing Seurat object for patient: ", patient_id)
    p <- preprocess_seurat(p, my_verbose, blueprint_encode)
    message("Completed preprocessing Seurat object for patient: ", patient_id)
    return(p)
  }, error = function(e) {
    stop("Error preprocessing Seurat object for patient: ",
         patient_id, ": ", e$message)
  })
})

message("Finished preprocessing seurat object(s).")

# ========================================================
# Step 2: Data Integration and Batch Correction
# ========================================================

# Check if Seurat objects are ready for integration
is_seurat_ready_integration(patient_seurat_list, patients)

# Identify integration anchors using only the common features
features_to_integrate <-
  SelectIntegrationFeatures(object.list = patient_seurat_list,
                            nfeatures = 4000)

# Prepare for integration
patient_seurat_list <-
  PrepSCTIntegration(object.list = patient_seurat_list,
                     anchor.features = features_to_integrate,
                     verbose = my_verbose)

# Run PCA on each Seurat object
patient_seurat_list <- lapply(patient_seurat_list, function(seurat_obj) {
  RunPCA(seurat_obj, features = features_to_integrate, verbose = my_verbose)
})

anchors <- FindIntegrationAnchors(object.list = patient_seurat_list,
                                  anchor.features = features_to_integrate,
                                  dims = 1:30, normalization.method = "SCT",
                                  reduction = "rpca", k.anchor = 20,
                                  verbose = my_verbose, reference = 1)

# Clean up memory by removing temporary objects
rm(patient_seurat_list, features_to_integrate)

#' @todo Ask doctor about the warnings here
integrated_seurat <- IntegrateData(anchorset = anchors,
                                   normalization.method = "SCT", dims = 1:30,
                                   verbose = my_verbose)

# Clean up memory by removing anchors
rm(anchors)

integrated_seurat$type <- factor(integrated_seurat$type,
                                 levels = c("Early", "Late"))

# ========================================================
# Step 3: Clustering
# ========================================================

# Run PCA
integrated_seurat <-
  RunPCA(integrated_seurat,
         features = VariableFeatures(object = integrated_seurat))

# Run UMAP
integrated_seurat <-
  RunUMAP(integrated_seurat, dims = 1:50, verbose = my_verbose)

# Find neighbors
integrated_seurat <-
  FindNeighbors(integrated_seurat, dims = 1:50, verbose = my_verbose)

# Find neighbors and clusters
integrated_seurat <-
  FindClusters(integrated_seurat, resolution = seq(0.1, 1, by = 0.1),
               verbose = my_verbose, algorithm = 1)

# Find the best resolution using silhouette scores
silhouette_results <-
  find_best_resolution(integrated_seurat, seq(0.1, 1, by = 0.1),
                       pca_dims = 1:50)
best_resolution <- silhouette_results$best_resolution

# Plot silhouette scores
plot_silhouette_scores(seq(0.1, 1, by = 0.1), silhouette_results$mean_scores,
                       silhouette_results$sd_scores, plot_output_path)

# Set clusters based on the best resolution
integrated_seurat <- set_best_clusters(integrated_seurat, best_resolution)

# Rename clusters based on predefined labels
cluster_labels <- c("CD8 T-cell", "CD4 T-cell 1", "Plasma Cells", "Tumor.1",
                    "Tumor.2", "Tumor.3", "Tregs", "B-cells", "Myeloid",
                    "Endothelial", "Fibroblast", "Tumor.4", "CD4 T-cell 2",
                    "Misc")
integrated_seurat <-
  plot_umap_clusters(integrated_seurat, cluster_labels, plot_output_path)

# Find top genes per cluster
top_genes <- find_top_genes(integrated_seurat)

plot_gene_heatmap(
  GetAssayData(integrated_seurat, assay = "SCT", layer = "scale.data"),
  integrated_seurat$seurat_clusters,
  top_genes$gene,
  n_top_genes_per_cluster = 5,
  scaled = FALSE,
  plot_output_path = plot_output_path
)

# Refine labels based on blueprint labels and p-values
integrated_seurat <- filter_blueprint_labels(integrated_seurat)

# Plot UMAP with refined labels
plot_umap_with_labels(integrated_seurat, plot_output_path)

# Save the integrated data
saveRDS(integrated_seurat,
        file = file.path(base_output_path, "colorectal_integrated.rds"))

# ========================================================
# Step 4: Generating Metacell Matrices
# ========================================================

metacell_matrices <-
  generate_metacell_matrices(integrated_seurat, base_output_path, "metacell")

# Plot cluster frequencies by treatment
plot_cluster_freq_by_treatment(integrated_seurat, plot_output_path)

# ========================================================
# Step 5: Running ARACNe and VIPER Analysis
# ========================================================

if (length(metacell_matrices) > 0) {
  expression_files <-
    prep_and_save_expr_for_aracne(base_output_path, "_all_all.txt.tsv")
  aracne_output_base_dir <- file.path(base_output_path, "aracne_results")
  run_aracne_for_all(aracne_binary_path, expression_files, regulator_files,
                     aracne_output_base_dir, threads = 4, seed = 42)

  # Load expression matrix from Seurat object
  exp_mat <-
    GetAssayData(object = integrated_seurat, assay = "SCT", layer = "data")
  if (is.null(exp_mat) || ncol(exp_mat) == 0 || nrow(exp_mat) == 0) {
    stop(paste0("Expression matrix is empty or NULL.",
                " Check your Seurat object and data extraction steps."))
  }
  exp_mat <- as.matrix(exp_mat)

  # Process ARACNe output files to generate regulon objects
  regulon_list <-
    generate_regulon_objects(aracne_output_base_dir, exp_mat, base_output_path)

  # Run VIPER analysis on the regulon objects
  viper_results <- run_viper(exp_mat, regulon_list)
  viper_results_path <- file.path(base_output_path, "viper_results.rds")
  save_viper_results(viper_results, viper_results_path)

  message("VIPER analysis completed and results saved.")
} else {
  cat(paste0("No metacell matrices were generated.",
             "Skipping ARACNe and analysis.\n"))
}

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ========================================================
# Step 6: Re-clustering based on VIPER results
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
# Step 7: Plotting the frequency of each cluster by patient
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