source("classes/utils.r")
source("classes/converter.r")
source("classes/loader.r")
source("classes/preprocessor.r")
source("classes/integrator.r")
source("classes/clusterer.r")
source("classes/plotter.r")
source("classes/metacell_generator.r")
source("classes/runner.r")

####################### DEFINE YOUR LOCAL PATHS HERE ##########################

#' @todo define the path of the directory where each of the patient directories
#' are located. For example if you had patient directories P1, P2,..., Pn. They
#' would be located at base_data_path/Pi.
#'
#' @note if by any chance you have data that's not in the Read10X format
#' (e.g. .rds ensemble), you can still use this script. Just make sure that the
#' .rds files are in the format base_data_path/Pi.rds and that you do step 0!
#'
#' @note: The directory should not contain any other directories/files than the
#' patient directories/files!
base_data_path <- paste0("/path/to/data")

#' @todo define the path of the of the data directory within each patient
#' directory.
#' To be precise imagine the following:
#' You cd into a patient directory, where would you find the data?
#' For example, if you cd into P1, you would find the data in P1/data. In this
#' case, `patient_data_path` would be "data".
#'
#' @note: The each patient directory must have the same path within that leads
#' to the data directory.
patient_data_path <- paste0("analysis/cellranger-count-default/",
                            "cellranger_count_outs/filtered_feature_bc_matrix")

#' @todo define the path of the directory where you want all your output to go.
base_output_path <- "/path/to/output"

#' @todo define the path of your ARACNe3 binary executable on the machine you
#' are running this script on.
#'
#' @note: It must be ARAcNe3, not ARACNe2.
aracne_binary_path <- paste0("/path/to/ARACNe3/build/src/app/",
                             "ARACNe3_app_release")

#' @todo define paths to regulator files, currently supported only in .txt
#' format
regulator_dir_path <- paste0("/path/to/regulators/human_hugo")
regulator_files <- list(
  cotfs = file.path(regulator_dir_path, "cotfs-hugo.txt"),
  surface = file.path(regulator_dir_path, "surface-hugo.txt"),
  sig = file.path(regulator_dir_path, "sig-hugo.txt"),
  tfs = file.path(regulator_dir_path, "tfs-hugo.txt")
)

############################## DEFINE METADATA ################################

#' @todo define any metadata you want.
#
#' @note: Please make sure you include an `id` in your metadata for
#' each patient. The patient id should be the name of the topmost directory
#' containing that patient's data. For example, if the data for patient P1 is
#' located at base_data_path/P1, then the id for P1 should be "P1".
patients <- list(
  list(id = "CRC0008", type = "Late"),
  list(id = "CRC0026", type = "Late"),
  list(id = "CRC0080", type = "Early"),
  list(id = "CRC0081", type = "Early"),
  list(id = "CRC0084", type = "Late"),
  list(id = "JD001", type = "Early"),
  list(id = "JD002", type = "Early"),
  list(id = "JD004", type = "Late"),
  list(id = "JD005", type = "Early"),
  list(id = "JD006", type = "Late")
)

########################### DEFINE OTHER PREFERENCES ##########################

#' @todo set to TRUE to display verbose messages during the analysis
my_verbose <- FALSE

# =============================================================================
# Initialization (Don't add/remove anything here)
# =============================================================================

utils <- Utils$new()
paths <- utils$init(base_data_path, base_output_path,
                    aracne_binary_path, regulator_dir_path)

plot_output_path <- paths$plot_output_path
plot_viper_output_path <- paths$plot_viper_output_path
aracne_output_path <- paths$aracne_output_path
viper_output_path <- paths$viper_output_path

rm(Utils, utils)
# =============================================================================
# Step 0: Conversion for Non-Read10X Data
# =============================================================================

#' @todo in case that your data is not in the Read10X format, you can convert
#' the data to the required format here. For example, if you have .rds files
#' with ensemble IDs, you can convert them to gene names like this:
# converter <- Converter$new(base_data_path)
# converter$convert_ensembl_to_gene_names()

# rm(Converter, converter)
# =============================================================================
# Step 1: Load data (No action necessary)
# =============================================================================

loader <- Loader$new(patients, base_data_path, patient_data_path)
patient_seurat_list <- loader$load_data()

# User can perform intermediate steps here if needed, such as plotting

rm(Loader, loader)
# =============================================================================
# Step 2: Preprocess data (No action necessary)
# =============================================================================

preprocessor <- Preprocessor$new(patient_seurat_list, verbose = my_verbose)
patient_seurat_list <- preprocessor$preprocess_data()

# User can perform intermediate steps here if needed, such as plotting

rm(Preprocessor, preprocessor)
# =============================================================================
# Step 3: Data Integration and Batch Correction (No action necessary)
# =============================================================================

integrator <- Integrator$new(patient_seurat_list, verbose = my_verbose)
integrated_seurat <- integrator$integrate_data()

rm(Integrator, integrator)
# =============================================================================
# Step 4: Clustering
# =============================================================================

#' @todo Define cluster labels if you know them.
cluster_labels <- c("CD8 T-cell", "CD4 T-cell 1", "Plasma Cells", "Tumor.1",
                    "Tumor.2", "Tumor.3", "Tregs", "B-cells", "Myeloid",
                    "Endothelial", "Fibroblast", "Tumor.4", "CD4 T-cell 2",
                    "Misc")

# Perform clustering
clusterer <- Clusterer$new(integrated_seurat, verbose = my_verbose)
clusterer$run_clustering()

# Find the best resolution based on silhouette scores
silhouette_results <- clusterer$calc_silhouette_scores()
best_resolution <- silhouette_results$best_resolution

# Set clusters based on the best resolution
integrated_seurat <- clusterer$set_best_clusters(best_resolution)

# Find top genes per cluster
top_genes <- clusterer$find_top_genes(assay_name = "SCT")

# Plotting
plotter <- Plotter$new(integrated_seurat, plot_output_path, cluster_labels)

# Plot silhouette scores
plotter$plot_silhouette_scores(silhouette_results$mean_scores,
                               silhouette_results$sd_scores)

# Plot UMAP clusters
plotter$plot_umap_clusters()

# Plot gene heatmap
plotter$plot_gene_heatmap(top_genes$gene, assay = "SCT")

# Plot UMAP with refined labels
plotter$plot_umap_with_labels()

integrated_seurat <- plotter$seurat_obj

# Save the integrated data
saveRDS(integrated_seurat,
        file = file.path(base_output_path, "integrated_seurat.rds"))

rm(Clusterer, clusterer)
# =============================================================================
# Step 5: Generating Metacell Matrices
# =============================================================================

generator <- MetacellGenerator$new(integrated_seurat, base_output_path)
metacell_matrices <- generator$generate_metacell_matrices()

#' @todo define the plot title and the order of the `type` groups
plot_title <- "Cluster Frequency by Early vs Late"
group_levels <- c("Early", "Late")

# Plot cluster frequencies
plotter$plot_cluster_freq_by(plot_title, group_by = "type", plot_type = "dot",
                             group_levels = group_levels)
plotter$plot_cluster_freq_by(plot_title, group_by = "type", plot_type = "box",
                             group_levels = group_levels)

rm(MetacellGenerator, generator)
# =============================================================================
# Step 6: Running ARACNe and VIPER Analysis (No action necessary)
# =============================================================================

runner <- Runner$new(
  metacell_matrices, base_output_path, aracne_binary_path, aracne_output_path,
  regulator_files, integrated_seurat, viper_output_path, threads = 4, seed = 42
)

runner$run_aracne()
runner$run_viper()

rm(Runner, runner)
# =============================================================================
# Step 7: Re-clustering Based on VIPER Results
# =============================================================================

viper_results <-
  readRDS(file = file.path(viper_output_path, "viper_results.rds"))

#' @note: The integrated Seurat object is reloaded here. If you've already
#' retained it in your environment, you may skip this step.
integrated_seurat <-
  readRDS(file = file.path(base_output_path, "integrated_seurat.rds"))

# Attach VIPER results as a new assay in the Seurat object
integrated_seurat[["VIPER"]] <- CreateAssayObject(counts = viper_results)

# Set the default assay to VIPER
DefaultAssay(integrated_seurat) <- "VIPER"

viper_features <- rownames(integrated_seurat[["VIPER"]])
if (length(viper_features) < 2) {
  stop("Not enough features in the VIPER assay to run PCA.")
}
VariableFeatures(integrated_seurat, assay = "VIPER") <- viper_features

# ScaleData is deliberately NOT called on the VIPER assay.
#
# aREA already returns NES: a z-like statistic, comparable across regulators and
# cells, sign-interpretable, positive = active. ScaleData z-scores each protein
# ACROSS cells, which is a second normalisation of an already-normalised
# quantity. It forces every regulator to mean 0 and sd 1, so a protein that is
# genuinely active in most cells is flattened to look average, and a uniformly
# inactive one is inflated into apparent structure. Measured on a 4,804-protein
# x 203,516-cell object, the SD of per-protein means went 0.214 -> 0.000: the
# baseline-activity differences that make protein activity worth computing are
# exactly what gets removed.
#
# The heatmap and PCA both read scale.data, so it still has to be populated -
# with raw NES rather than a rescaling of it.
nes <- LayerData(integrated_seurat, assay = "VIPER", layer = "counts")
integrated_seurat <- SetAssayData(integrated_seurat, assay = "VIPER",
                                  layer = "data", new.data = nes)
integrated_seurat <- SetAssayData(integrated_seurat, assay = "VIPER",
                                  layer = "scale.data",
                                  new.data = as.matrix(nes))
rm(nes)

# Perform PCA on the VIPER assay
integrated_seurat <- RunPCA(integrated_seurat,
                            assay = "VIPER",
                            features = viper_features,
                            verbose = my_verbose)

# Re-cluster the data based on VIPER results
clusterer_viper <- Clusterer$new(integrated_seurat, verbose = my_verbose)
clusterer_viper$run_clustering()

silhouette_results_viper <- clusterer_viper$calc_silhouette_scores()
best_resolution_viper <- silhouette_results_viper$best_resolution

integrated_seurat <- clusterer_viper$set_best_clusters(best_resolution_viper)

top_regulators_viper <- clusterer_viper$find_top_regulators()
write.csv(top_regulators_viper$all,
          file.path(viper_output_path, "viper_cluster_markers.csv"),
          row.names = FALSE)

saveRDS(integrated_seurat,
        file = file.path(base_output_path,
                         "integrated_seurat_viper_reclustered.rds"))

#' @note: Congrats! You've successfully re-clustered your data based on VIPER
#' results. You can now proceed with plotting. You wonder how to do that?
#' Here's a template for you:

# You can define cluster labels specifically for your VIPER-based
# clustering here. For example:
cluster_labels_viper <- c("VIPER_C1", "VIPER_C2", "VIPER_C3", ...)

# For plotting you can use the same Plotter class as before, but with the new
# cluster labels and the re-clustered Seurat object.
plotter_viper <- Plotter$new(integrated_seurat,
                             plot_viper_output_path,
                             cluster_labels_viper)

rm(Clusterer, clusterer_viper, Plotter, plotter_viper)