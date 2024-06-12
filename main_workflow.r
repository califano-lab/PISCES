source("classes/loader.r")
source("classes/preprocessor.r")
source("classes/integrator.r")
source("classes/clusterer.r")
source("classes/plotter.r")
source("classes/metacell_generator.r")
source("classes/runner.r")
source("classes/utils.r")

################## DEFINE YOUR LOCAL PATHS HERE ##################

#' @todo define the path of the directory where each of the patient directories
#' are located. For example if you had patient directories P1, P2,..., Pn. They
#' would be located at base_data_path/Pi.
#'
#' @note: The directory should not contain any other directories than the
#' patient directories!
base_data_path <- paste0("/Users/apple/Documents/Research/aleks-lab/",
                         "data/colorectal-data")

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
base_output_path <- "/Users/apple/Desktop/test-output"

#' @todo define the path of your ARACNe3 binary executable on the machine you
#' are running this script on.
#'
#' @note: It must be ARAcNe3, not ARACNe2.
aracne_binary_path <- paste0("/Users/apple/Documents/Research/aleks-lab/",
                             "repos/ARACNe3/build/src/app/",
                             "ARACNe3_app_release")

#' @todo define paths to regulator files, currently supported only in .txt
#' format
regulator_dir_path <- paste0("/Users/apple/Documents/Research/aleks-lab/",
                             "data/regulators/human_hugo")
regulator_files <- list(
  cotfs = file.path(regulator_dir_path, "cotfs-hugo.txt"),
  surface = file.path(regulator_dir_path, "surface-hugo.txt"),
  sig = file.path(regulator_dir_path, "sig-hugo.txt"),
  tfs = file.path(regulator_dir_path, "tfs-hugo.txt")
)

################## DEFINE METADATA ##############################

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

################## DEFINE OTHER PREFERENCES #####################

#' @todo set to TRUE to display verbose messages during the analysis
my_verbose <- FALSE

# =============================================================================
# Step 0: Initialization (Don't add/remove anything here)
# =============================================================================

utils <- Utils$new()
paths <- utils$init(base_output_path, plot_output_path, aracne_output_path,
                    viper_output_path, base_data_path, aracne_binary_path,
                    regulator_dir_path)

plot_output_path <- paths$plot_output_path
aracne_output_path <- paths$aracne_output_path
viper_output_path <- paths$viper_output_path

rm(Utils, utils)
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
top_genes <- clusterer$find_top_genes()

# Plotting
plotter <- Plotter$new(integrated_seurat, plot_output_path, cluster_labels)

# Plot silhouette scores
plotter$plot_silhouette_scores(silhouette_results$mean_scores,
                               silhouette_results$sd_scores)

# Plot UMAP clusters
plotter$plot_umap_clusters()

# Plot gene heatmap
plotter$plot_gene_heatmap(top_genes$gene)

# Plot UMAP with refined labels
plotter$plot_umap_with_labels()

integrated_seurat <- plotter$seurat_obj

# Save the integrated data
saveRDS(integrated_seurat,
        file = file.path(base_output_path, "colorectal_integrated.rds"))

rm(Clusterer, clusterer)
# =============================================================================
# Step 5: Generating Metacell Matrices
# =============================================================================

generator <- MetacellGenerator$new(integrated_seurat, base_output_path)
metacell_matrices <- generator$generate_metacell_matrices()

#' @todo define the column names and plot title for the cluster frequency plot
col_names <-
  c("Early_p1", "Early_p2", "Early_p3", "Late_p1", "Late_p2", "Late_p3")
plot_title <- "Cluster Frequency by Early vs Late"

# Plot cluster frequencies by treatment
plotter$plot_cluster_freq_by_treatment(col_names, plot_title)

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

rm(Runner, runner, Plotter, plotter)
# =============================================================================
# Step 7, ...: Re-clustering based on VIPER results, ... (SOON)
# =============================================================================