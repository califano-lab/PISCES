source("loader.r")
source("preprocessor.r")
source("integrator.r")
source("clusterer.r")
source("plotter.r")
source("metacell_generator.r")
source("runner.r")
source("aracne.r")
source("viper.r")
source("utils.r")

library(dplyr)
library(ggplot2)

################## DEFINE YOUR LOCAL PATHS HERE ##################

# Define the path of the directory where each of the patient directories
# are located. For example if you had patient directories P1, P2,..., Pn. They
# would be located at base_data_path/Pi.
#
# Very important: The directory should not contain any other directories than
# the patient directories!
base_data_path <- "/Users/apple/Desktop/colorectal-data"

# Define the path of the of the data directory within each patient directory.
# To be precise imagine the following:
# You cd into a patient directory, where would you find the data?
# For example, if you cd into P1, you would find the data in P1/data. In this
# case, data_dir_path would be "data".
#
# Very important: The each patient directory must have the same path within
# that leads to the data directory.
patient_data_path <- paste0("analysis/cellranger-count-default/",
                            "cellranger_count_outs/filtered_feature_bc_matrix")

# Define the path of the directory where you want all your output to go.
base_output_path <- "/Users/apple/Desktop/test-output"

# Define the path of your ARACNe3 binary executable on the machine you are
# running this script on.
aracne_binary_path <- paste0("/Users/apple/Documents/Research/aleks-lab/",
                             "repos/ARACNe3/build/src/app/",
                             "ARACNe3_app_release")

# Define paths to regulator files, currently supported only in .txt format
regulator_dir_path <- "/Users/apple/Desktop/output/human_hugo"
regulator_files <- list(
  cotfs = file.path(regulator_dir_path, "cotfs-hugo.txt"),
  surface = file.path(regulator_dir_path, "surface-hugo.txt"),
  sig = file.path(regulator_dir_path, "sig-hugo.txt"),
  tfs = file.path(regulator_dir_path, "tfs-hugo.txt")
)

plot_output_path <- file.path(base_output_path, "plots")
aracne_output_path <- file.path(base_output_path, "aracne_results")
viper_output_path <- file.path(base_output_path, "viper_results")

utils <- Utils$new()
# Create necessary directories
utils$create_directories(c(
  base_output_path,
  plot_output_path,
  aracne_output_path,
  viper_output_path
))

utils$do_directories_exist(
  c(base_data_path, base_output_path, aracne_binary_path, regulator_dir_path)
)
rm(utils)

#################################################################

################## DEFINE METADATA ##############################

# Define any metadata you want.
#
# Very important: Please make sure you include an `id` in your metadata for
# each patient. The patient id should be the name of the topmost directory
# containing that patient's data. For example, if the data for patient P1 is
# located at base_data_path/P1, then the id for P1 should be "P1".
patients <- list(
  list(id = "JD001", type = "Early"),
  list(id = "JD002", type = "Early"),
  list(id = "JD003", type = "Late"),
  list(id = "JD004", type = "Late"),
  list(id = "JD005", type = "Early"),
  list(id = "JD006", type = "Late")
)

################## DEFINE OTHER PREFERENCES #####################

# Set to TRUE to display verbose messages during the analysis
my_verbose <- FALSE

#################################################################

# ========================================================
# Step 1: Load data
# ========================================================

loader <- Loader$new(patients, base_data_path, patient_data_path)
patient_seurat_list <- loader$load_data()
rm(Loader, loader)

# User can perform intermediate steps here if needed, such as plotting

# ========================================================
# Step 2: Preprocess data
# ========================================================

preprocessor <- Preprocessor$new(patient_seurat_list, verbose = my_verbose)
patient_seurat_list <- preprocessor$preprocess_data()
rm(Preprocessor, preprocessor)

# User can perform intermediate steps here if needed, such as plotting

# ========================================================
# Step 3: Data Integration and Batch Correction
# ========================================================

integrator <- Integrator$new(patient_seurat_list, verbose = my_verbose)
integrated_seurat <- integrator$integrate_data()

# integrated_seurat$type <- factor(integrated_seurat$type,
#                                  levels = c("Early", "Late"))

rm(Integrator, integrator)
# ========================================================
# Step 4: Clustering
# ========================================================

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

# $$$ Sort of up to here $$$

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
# ========================================================
# Step 5: Generating Metacell Matrices
# ========================================================

generator <- MetacellGenerator$new(integrated_seurat, base_output_path)
metacell_matrices <- generator$generate_metacell_matrices()

col_names <-
  c("Early_p1", "Early_p2", "Early_p3", "Late_p1", "Late_p2", "Late_p3")
plot_title <- "Cluster Frequency by Early vs Late"

# Plot cluster frequencies by treatment
plotter$plot_cluster_freq_by_treatment(col_names, plot_title)

rm(MetacellGenerator, generator)
# ========================================================
# Step 6: Running ARACNe and VIPER Analysis
# ========================================================

runner <- Runner$new(
  metacell_matrices, base_output_path, aracne_binary_path, aracne_output_path,
  regulator_files, integrated_seurat, viper_output_path, threads = 4, seed = 42
)

runner$run_aracne()
runner$run_viper()

rm(Runner, runner)