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

# Create the directory if it does not exist
if (!dir.exists(plot_output_path)) {
  dir.create(plot_output_path, recursive = TRUE)
}
#################################################################

################## DEFINE OTHER PREFERENCES #####################

my_verbose <- FALSE

#################################################################

# ========================================================
# Step 1: Load and process data
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

# Load data for each patient
patient_data_list <- lapply(patients, function(patient) {
  tryCatch({
    message("Loading patient: ", patient$id)
    seurat_obj <- load_patient_data(patient, base_path, analysis_prefix,
                                    count_default_suffix, output_folder_suffix,
                                    feature_matrix_dir)
    message("Completed loading for patient: ", patient$id)
    return(seurat_obj)
  }, error = function(e) {
    message("Error loading patient: ", patient$id, ": ", e$message)
    return(NULL)
  })
})

message("Finished loading patient data.")

# Filter out any NULL values that resulted from errors
patient_data_list <- Filter(Negate(is.null), patient_data_list)

blueprint_encode <- BlueprintEncodeData()

patient_data_list <- lapply(patient_data_list, function(p) {
  patient_id <- unique(p$patient)
  tryCatch({
    message("Processing Seurat object for patient: ", patient_id)
    p <- process_patient_data(p, my_verbose, blueprint_encode)
    message("Completed processing Seurat object for patient: ", patient_id)
    return(p)
  }, error = function(e) {
    message("Error processing Seurat object for patient: ",
            patient_id, ": ", e$message)
    return(NULL)
  })
})


message("Finished processing patient data.")

# Filter out any NULL values that resulted from errors
patient_data_list <- Filter(Negate(is.null), patient_data_list)

# ========================================================
# Step 2: Data Integration and Batch Correction
# ========================================================

# Check if Seurat objects are ready for integration
invisible(lapply(names(patient_data_list), function(patient_id) {
  seurat_object <- patient_data_list[[patient_id]]
  if ("RNA" %in% names(seurat_object@assays) &&
        ncol(GetAssayData(seurat_object, assay = "RNA",
                          layer = "counts")) > 0) {
    message(paste("Seurat object for patient", patient_id,
                  "is ready for integration."))
  } else {
    stop(paste("Seurat object for patient", patient_id,
               "is not ready for integration."))
  }
}))

# Identify integration anchors using only the common features
features_to_integrate <-
  SelectIntegrationFeatures(object.list = patient_data_list, nfeatures = 4000)

# Prepare for integration
patient_data_list <-
  PrepSCTIntegration(object.list = patient_data_list,
                     anchor.features = features_to_integrate,
                     verbose = my_verbose)

patient_data_list <- lapply(patient_data_list, FUN = RunPCA,
                            features = features_to_integrate)

anchors <- FindIntegrationAnchors(object.list = patient_data_list,
                                  anchor.features = features_to_integrate,
                                  dims = 1:30, normalization.method = "SCT",
                                  reduction = "rpca", k.anchor = 20,
                                  verbose = my_verbose, reference = 1)

# Clean up memory by removing temporary objects
rm(patient_data_list, features_to_integrate)

integrated_data <- IntegrateData(anchorset = anchors,
                                 normalization.method = "SCT", dims = 1:30,
                                 verbose = my_verbose)

# Clean up memory by removing anchors
rm(anchors)

integrated_data$type <- factor(integrated_data$type,
                               levels = c("Early", "Late"))

# ========================================================
# Step 3: Clustering and Identifying Regulatory Networks
# ========================================================

# Running PCA on the integrated data to enable visualization and
# further analysis
integrated_data <- RunPCA(integrated_data, verbose = TRUE)

# Generate a UMAP reduction for visualization
integrated_data <- RunUMAP(integrated_data, reduction = "pca", dims = 1:20)

p <- DimPlot(integrated_data, reduction = "umap", group.by = "patient",
             label = TRUE) +
  ggtitle("UMAP Visualization of Integrated Single-cell Data") +
  scale_color_viridis_d() +
  theme(legend.position = "right")

# Define the path and filename for the UMAP plot
umap_plot_path <- file.path(plot_output_path, "umap_integration_results.png")

# Save the UMAP plot
ggsave(umap_plot_path, plot = p, width = 10, height = 8)
cat("UMAP plot saved to:", umap_plot_path, "\n")



# Find neighbors and clusters
integrated_data <- FindNeighbors(integrated_data, dims = 1:20)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)

# UMAP plot of clusters
p <- DimPlot(integrated_data, reduction = "umap",
             group.by = "seurat_clusters") +
  ggtitle("UMAP Clustering Results")

# Define path for saving the UMAP plot
umap_cluster_plot_path <-
  file.path(plot_output_path, "umap_clustering_results.png")
ggsave(umap_cluster_plot_path, plot = p, width = 10, height = 8)

# Verify the RNA assay and set it as the default
if ("RNA" %in% names(integrated_data@assays)) {
  DefaultAssay(integrated_data) <- "RNA"

  # Find top genes per cluster with extremely lenient thresholds
  top_genes <- FindAllMarkers(integrated_data, only.pos = TRUE, min.pct = 0.1,
                              logfc.threshold = 0.5, test.use = "wilcox")

  if (nrow(top_genes) > 0) {
    # Group and select top 5 genes per cluster by log fold change
    top_genes <- top_genes %>%
      dplyr::group_by(cluster) %>%
      dplyr::top_n(n = 5, wt = avg_log2FC)

    if (nrow(top_genes) > 0) {
      gene_list <- top_genes$gene
      data_matrix <-
        GetAssayData(integrated_data,
                     layer = "data")[gene_list, , drop = FALSE]

      if (!is.null(data_matrix) && ncol(data_matrix) > 0 &&
            length(gene_list) == nrow(data_matrix)) {
        # Generate heatmap
        heatmap_plot <-
          gene_heatmap_plot(data_matrix,
                            integrated_data@meta.data$seurat_clusters,
                            genes = gene_list, genes_by_cluster = TRUE,
                            n_top_genes_per_cluster = 5, scaled = TRUE)

        # Save the heatmap
        heatmap_plot_path <-
          file.path(plot_output_path, "gene_expression_heatmap.png")
        ggsave(heatmap_plot_path, plot = heatmap_plot, width = 10, height = 8)
      } else {
        message(paste0("Heatmap data matrix is not valid for plotting.",
                       " Check gene list and data matrix dimensions."))
      }
    } else {
      message(paste0("Not enough significant markers",
                     " found after grouping by cluster."))
    }
  } else {
    message(paste0("No significant markers found across any clusters.",
                   " Consider adjusting the thresholds or revising",
                   " the clustering approach."))
  }
} else {
  message("RNA assay not found in the dataset.")
}

# Ensure the SCT assay is set as default if it's being used
DefaultAssay(integrated_data) <- "SCT"

# Try to access the normalized data from the SCT assay
if ("SCT" %in% names(integrated_data@assays)) {
  counts_matrix <-
    GetAssayData(object = integrated_data, assay = "SCT", layer = "data")
  if (is.null(counts_matrix) || ncol(counts_matrix) == 0 ||
        nrow(counts_matrix) == 0) {
    stop(paste0("Normalized data matrix is empty or NULL.",
                "Check your Seurat object and data extraction steps."))
  } else {
    print(paste("Normalized data matrix dimensions:", nrow(counts_matrix),
                "genes X", ncol(counts_matrix), "cells"))
  }
} else {
  stop("SCT assay not found in the integrated Seurat object.")
}

# Call the function with the correctly obtained counts matrix
metacell_matrices <- make_cmfa(dat_mat = counts_matrix,
                               clustering = integrated_data@active.ident,
                               out_dir = base_output_path,
                               out_name = "metacell")

if (length(metacell_matrices) > 0) {
  # Prepare and save the expression data for each cluster
  expression_files <-
    prep_and_save_expr_for_aracne(base_output_path, "_all_all.txt.tsv")

  # Get regulators from a file
  regulators_file_path <- paste0(base_output_path, "/regulators.txt")
  regulators <- readLines(paste0(base_output_path, "/regulators_list.txt"))
  writeLines(regulators, con = regulators_file_path)

  # Running ARACNe for each cluster file
  aracne_output_dir <- paste0(base_output_path, "/aracne_results")
  lapply(expression_files, function(exp_file) {
    run_aracne(aracne_binary_path, exp_file, regulators_file_path,
               aracne_output_dir, threads = 4, seed = 42)
  })

  # Integrating VIPER analysis
  # Process for ARACNe and VIPER
  aracne_output_files <- list.files(aracne_output_dir,
                                    pattern = "consolidated-net_.*\\.tsv$",
                                    full.names = TRUE)

  viper_results <- lapply(aracne_output_files, function(aracne_file) {
    cat("Processing ARACNe output file:", aracne_file, "\n")

    # Load ARACNe output file without headers
    aracne_data <- read.table(aracne_file, header = FALSE, sep = "\t",
                              check.names = FALSE, stringsAsFactors = FALSE,
                              skip = 1)

    #Only include the first three columns
    aracne_data <- aracne_data[, 1:3]

    # Define column names manually
    colnames(aracne_data) <- c("regulator", "target", "mi")

    # Convert the 'mi' column to numeric, checking for invalid data
    aracne_data$mi <- as.numeric(aracne_data$mi)

    # Load expression matrix from Seurat object
    exp_mat <-
      GetAssayData(object = integrated_data, assay = "SCT", slot = "data")
    if (is.null(exp_mat) || ncol(exp_mat) == 0 || nrow(exp_mat) == 0) {
      stop(paste0("Expression matrix is empty or NULL.",
                  " Check your Seurat object and data extraction steps."))
    }
    print(paste("Expression matrix dimensions: Genes =", nrow(exp_mat),
                "Samples =", ncol(exp_mat)))

    ## Convert exp_mat to a matrix
    exp_mat <- as.matrix(exp_mat)

    # Convert and analyze
    regulon_object <- convert_to_regulon(aracne_data, exp_mat)
    viper_scores <- tryCatch({
      viper::viper(exp_mat, regulon_object)
    }, error = function(e) {
      cat("Error in VIPER analysis:", e$message, "\n")
      NULL
    })

    return(viper_scores)
  })

  # Debugging: check contents of viper_results
  if (length(viper_results) == 0 || any(sapply(viper_results, is.null))) {
    stop("VIPER results are empty or not properly formed.")
  }

  # Save VIPER results
  viper_results_path <- paste0(base_output_path, "/viper_results.rds")
  saveRDS(viper_results, file = viper_results_path)
} else {
  cat(paste0("No metacell matrices were generated.",
             "Skipping ARACNe and VIPER analysis.\n"))
}

# ========================================================
# Step 4: Re-clustering based on VIPER results
# ========================================================

# Verify that viper_results is a list and correctly formatted
if (!is.list(viper_results)) {
  stop("Error: VIPER results must be a list format.")
}

# Combine the list elements of viper_results into a matrix
viper_data <- do.call(cbind, viper_results)
colnames(viper_data) <- colnames(integrated_data)

# Output dimensions of viper_data to verify alignment with expected structure
cat("Dimensions of combined VIPER data matrix:")
print(dim(viper_data))

# The # of columns in viper_data matches the # of cells in integrated_data?
if (ncol(viper_data) != ncol(integrated_data)) {
  stop("Error: Mismatch in # of cells between VIPER and integrated data.",
       " Expected", ncol(integrated_data), "cells, but got",
       ncol(viper_data), "in VIPER results.")
}

# Validate that the column names in viper_data and integrated_data match
if (!all(colnames(viper_data) == colnames(integrated_data))) {
  stop(paste0("Error: Column names of the VIPER data do not",
              " match those of the Seurat object."))
}

# Add VIPER scores as a new assay in the Seurat object
integrated_data[["VIPER_scores"]] <- CreateAssayObject(viper_data)
DefaultAssay(integrated_data) <- "VIPER_scores"

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
integrated_data[["VIPER_scores"]] <- CreateAssayObject(non_constant_features)
DefaultAssay(integrated_data) <- "VIPER_scores"

# Perform hierarchical clustering on transposed data
dist_matrix <- dist(t(non_constant_features))
hc <- hclust(dist_matrix)
clusters <- cutree(hc, k = 5)

# Add cluster assignments to Seurat object
integrated_data <- AddMetaData(integrated_data, metadata = clusters,
                               col.name = "seurat_clusters")

# Visualize the re-clustering results using UMAP
p <- DimPlot(integrated_data, reduction = "umap",
             group.by = "seurat_clusters") +
  ggtitle("Re-clustering Based on VIPER Scores")

# Define the path and filename for saving the UMAP plot of re-clustered data
reclustered_umap_plot_path <-
  file.path(plot_output_path, "reclustered_umap_results.png")
ggsave(reclustered_umap_plot_path, plot = p, width = 10, height = 8)

# ========================================================
# Step 5: Plotting the frequency of each cluster by patient
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
p1 <- plot_cluster_frequency(integrated_data,
                             cluster_label = "seurat_clusters",
                             patient_label = "patient")

# Define the path and filename for saving the cluster frequency plot
cluster_frequency_plot_path <-
  file.path(plot_output_path, "cluster_frequency_by_patient.png")
ggsave(cluster_frequency_plot_path, plot = p1, width = 10, height = 8)

# Save the final integrated and annotated Seurat object
saveRDS(integrated_data,
        file = file.path(base_output_path,
                         "final_integrated_seurat_object.rds"))

# Identify top genes per cluster if not predefined
top_genes_per_cluster <- FindAllMarkers(integrated_data, only.pos = TRUE,
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
  gene_data <- GetAssayData(integrated_data, layer = "data")
  if (any(!top_genes$gene %in% rownames(gene_data))) {
    cat("Some genes in top_genes not found in the data matrix.\n")
  } else {
    p2 <- gene_heatmap_plot(gene_data,
                            integrated_data$seurat_clusters,
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