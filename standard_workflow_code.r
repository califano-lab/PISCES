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
library(celldex)
library(reshape2)
library(readr)

#' Create a large SingleR object with improved control over tuning
#'
#' This function creates a SingleR object for large datasets by processing in
#' chunks, addressing the issue in the original SingleR where fine tuning was
#' always enabled.
#'
#' @param counts                Expression count matrix with genes in rows and
#'                              cells in columns.
#' @param annot                 Optional annotation vector for cells.
#' @param project.name          Name of the project for file naming.
#' @param xy                    Coordinates for cells, typically used for
#'                              visualization.
#' @param clusters              Cluster assignments for cells.
#' @param N                     Number of cells to process in each chunk.
#' @param min.genes             Minimum number of genes for filtering cells.
#' @param technology            Single-cell technology used (e.g., "10X").
#' @param species               Species name (e.g., "Human").
#' @param citation              Citation or reference for the dataset.
#' @param ref.list              Reference list for SingleR classification.
#' @param normalize.gene.length Logical indicating whether to normalize by gene
#'                              length.
#' @param variable.genes        Method for selecting variable genes, default is
#'                              "de" (differentially expressed).
#' @param fine.tune             Logical indicating whether to fine-tune SingleR
#'                              results.
#' @param reduce.file.size      Logical indicating whether to reduce file size
#'                              by omitting some SingleR slots.
#' @param do.signatures         Logical indicating whether to compute signature
#'                              scores.
#' @param do.main.types         Logical indicating whether to classify main
#'                              cell types.
#' @param temp.dir              Directory to store temporary SingleR objects.
#' @param numCores              Number of cores to use for parallel processing.
#'
#' @return                      A combined SingleR object created from chunks.
create_big_single_r_object <- function(counts, annot = NULL, project_name,
                                       xy, clusters, n = 10000,
                                       min_genes = 200, technology = "10X",
                                       species = "Human", citation = "",
                                       ref_list = list(),
                                       normalize_gene_length = FALSE,
                                       variable_genes = "de",
                                       fine_tune = TRUE,
                                       reduce_file_size = TRUE,
                                       do_signatures = FALSE,
                                       do_main_types = TRUE,
                                       temp_dir = getwd(), num_cores = 1) {
  setup_temp_dir(temp_dir, project_name)
  cell_indices <- split_into_chunks(ncol(counts), n)

  for (A in cell_indices) {
    process_chunk(A, counts, annot, project_name, min_genes, technology,
                  species, citation, do_signatures, num_cores, fine_tune,
                  temp_dir)
  }

  singler_objects <- load_singler_objects(temp_dir, project_name)
  combine_singler_objects(singler_objects, colnames(counts), clusters, xy)
}

#' Setup Temporary Directory for SingleR Objects
#'
#' Creates a temporary directory within the specified path to store SingleR
#' object chunks.
#'
#' @param temp_dir     The base directory to create a temporary directory in.
#' @param project_name The name of the project, used in naming the temporary
#'                     directory.
setup_temp_dir <- function(temp_dir, project_name) {
  dir.create(file.path(temp_dir, "singler_temp"), showWarnings = FALSE)
}

#' Split Total Cells into Chunks
#'
#' Divides the total number of cells into smaller chunks for processing.
#'
#' @param total_cells Total number of cells in the dataset.
#' @param chunk_size  Desired number of cells in each chunk.
#'
#' @return            A list where each element contains cell indices for a
#'                    chunk.
split_into_chunks <- function(total_cells, chunk_size) {
  split(seq_len(total_cells), ceiling(seq_len(total_cells) / chunk_size))
}

#' Process a Chunk of Cells
#'
#' Processes a subset of cells to create a SingleR object for that chunk.
#'
#' @param cell_indices  Indices of cells in the chunk.
#' @param counts        Expression counts matrix.
#' @param annot         Cell annotations.
#' @param project_name  Name of the project.
#' @param min_genes     Minimum number of genes for inclusion.
#' @param technology    Single-cell sequencing technology used.
#' @param species       Species of the samples.
#' @param citation      Citation for the dataset.
#' @param do_signatures Whether to compute signature scores.
#' @param num_cores     Number of cores to use for computation.
#' @param fine_tune     Whether to fine-tune the SingleR results.
#' @param temp_dir      Temporary directory for storing SingleR objects.
process_chunk <- function(cell_indices, counts, annot, project_name, min_genes,
                          technology, species, citation, do_signatures,
                          num_cores, fine_tune, temp_dir) {
  singler <- SingleR::CreateSinglerObject(counts[, cell_indices],
                                          annot = annot[cell_indices],
                                          project_name = project_name,
                                          min_genes = min_genes,
                                          technology = technology,
                                          species = species,
                                          citation = citation,
                                          do_signatures = do_signatures,
                                          clusters = NULL,
                                          num_cores = num_cores,
                                          fine_tune = fine_tune)
  save(singler, file = file.path(temp_dir, "singler_temp",
                                 paste0(project_name, ".",
                                        cell_indices[1], ".RData")))
}

#' Load SingleR Objects from Files
#'
#' Loads SingleR objects saved in temporary files into a list.
#'
#' @param temp_dir     Directory where SingleR temporary files are stored.
#' @param project_name Name of the project, used to identify relevant files.
#'
#' @return             A list of SingleR objects.
load_singler_objects <- function(temp_dir, project_name) {
  singler_files <- list.files(file.path(temp_dir, "singler_temp"),
                              pattern = "RData", full.names = TRUE)

  lapply(singler_files, function(f) {
    # Load the .RData file
    loaded_names <- load(f)

    # Assume the object of interest is the last one loaded
    last_loaded_name <- tail(loaded_names, n = 1)

    # Use get() to retrieve the last loaded object by name
    get(last_loaded_name, envir = .GlobalEnv)
  })
}

#' Combine Singler objects into one
#'
#' @param singler_objects A list of Singler objects to combine.
#' @param cell_order      The order of cells to be maintained in the combined
#'                        object.
#' @param clusters        Cluster assignments for the cells.
#' @param xy              Coordinates for the cells, typically used for
#'                        visualization.
#'
#' @return                A combined Singler object.
combine_singler_objects <- function(singler_objects, cell_order,
                                    clusters, xy) {
  SingleR::SingleR.Combine(singler_objects, order = cell_order,
                           clusters = clusters, xy = xy)
}

#' Identify and Merge Unique Master Regulators for Each Cell
#'
#' This function computes the top master regulators (MRs) based on protein
#' activity for each cell in the input dataset. It then merges these lists and
#' returns a unique set of MRs across all cells.
#'
#' @param dat_mat A matrix with proteins as rows and samples (cells) as
#'                columns.
#' @param num_mrs Number of top MRs to identify in each cell.
#'
#' @return        A vector of unique master regulators identified across all
#'                cells.
cbcmrs <- function(dat_mat, num_mrs = 25) {
  # Identify MRs for each cell and merge into a unique list
  unique_mrs <- identify_and_merge_mrs(dat_mat, num_mrs)
  return(unique_mrs)
}

#' Identify and Merge MRs from Protein Activity Matrix
#'
#' This helper function applies over the columns of the protein activity matrix
#' to identify top MRs for each cell and then merges these lists, ensuring
#' uniqueness.
#'
#' @param dat_mat A matrix with proteins as rows and samples (cells) as
#'                columns.
#' @param num_mrs Number of top MRs to identify in each cell.
#'
#' @return        A vector of unique master regulators identified across all
#'                cells.
identify_and_merge_mrs <- function(dat_mat, num_mrs) {
  cbc_mrs <- apply(dat_mat, 2, function(x) {
    names(sort(x, decreasing = TRUE))[1:num_mrs]
  })
  unique(unlist(cbc_mrs))
}

#' Process ARACNe Results for VIPER Analysis
#'
#' Converts ARACNe output into a regulon object suitable for VIPER analysis,
#' including an optional pruning step to refine the regulon.
#'
#' @param a_file   Path to the ARACNe final network file in TSV format.
#' @param exp_mat  Expression matrix used to generate the ARACNe network
#'                 (genes x samples).
#' @param out_dir  Directory where the processed regulon objects will be saved.
#' @param out_name Prefix for the saved regulon files.
reg_process <- function(a_file, exp_mat, out_dir, out_name = "") {
  require(viper)

  # Convert ARACNe output to regulon object
  processed_reg <- convert_to_regulon(a_file, exp_mat)

  # Save the unpruned regulon object
  save_regulon(processed_reg, out_dir, paste0(out_name, "unpruned.rds"))

  # Prune the regulon to refine it
  pruned_reg <- prune_regulon(processed_reg)

  # Save the pruned regulon object
  save_regulon(pruned_reg, out_dir, paste0(out_name, "pruned.rds"))
}

#' Convert ARACNe Output to Regulon Object
#'
#' @param a_file  Path to the ARACNe network file.
#' @param exp_mat Expression matrix associated with the ARACNe network.
#'
#' @return        A regulon object suitable for VIPER analysis.
convert_to_regulon <- function(aracne_data, exp_mat) {
  if (is.null(aracne_data) || nrow(aracne_data) == 0) {
    stop("ARACNe data is empty or not available.")
  }

  if (is.null(dim(exp_mat))) {
    stop("Expression matrix is not correctly formatted or is NULL.")
  }

  # Create a temporary file to store processed ARACNe data
  temp_file <- tempfile()
  write.table(aracne_data, temp_file, sep = "\t", row.names = FALSE,
              col.names = FALSE, quote = FALSE)

  tryCatch({
    regulon_object <- aracne2regulon(afile = temp_file, eset = exp_mat,
                                     format = "3col", verbose = TRUE)
  }, error = function(e) {
    cat("Error during regulon conversion: ", e$message, "\n")
    stop("Failed to convert ARACNe output to regulon object: ", e$message)
  })

  unlink(temp_file)

  if (is.null(regulon_object) || length(regulon_object) == 0) {
    stop("Regulon object is NULL or empty.")
  }

  return(regulon_object)
}

#' Save Regulon Object to File
#'
#' @param regulon   The regulon object to be saved.
#' @param out_dir   The directory for saving the regulon file.
#' @param file_name The name of the file to save the regulon object in.
save_regulon <- function(regulon, out_dir, file_name) {
  saveRDS(regulon, file = file.path(out_dir, file_name))
}

#' Prune Regulon Object
#'
#' Applies pruning to a regulon object to refine its content.
#'
#' @param regulon The regulon object to be pruned.
#'
#' @return        A pruned regulon object.
prune_regulon <- function(regulon) {
  viper::pruneRegulon(regulon, 50, adaptive = FALSE, eliminate = TRUE)
}

#' Combine P-Values Using Fisher's Method
#'
#' This function combines multiple p-values into a single p-value using
#' Fisher's method. It's useful for meta-analyses where p-values from different
#' studies or tests need to be integrated.
#'
#' @param p_values A numeric vector of p-values to combine.
#'
#' @return         A single combined p-value.
combine_p_values <- function(p_values) {
  # Calculate the combined chi-squared statistic
  chi_squared_stat <- -2 * sum(log(p_values))

  # Determine the degrees of freedom for the test
  degrees_of_freedom <- 2 * length(p_values)

  # Calculate and return the combined p-value
  combined_p_value <- pchisq(chi_squared_stat, df = degrees_of_freedom,
                             lower.tail = FALSE)

  return(combined_p_value)
}

#' Combine Z-Scores Using Stouffer's Method
#'
#' This function combines multiple z-scores into a single z-score using
#' Stouffer's method. Optionally, it can weight the z-scores before combining
#' them.
#'
#' @param z_scores A numeric vector of z-scores to combine.
#' @param weights  Optional numeric vector of weights for each z-score;
#'                 defaults to equal weighting.
#'
#' @return         A single combined z-score.
combine_z_scores <- function(z_scores, weights = NULL) {
  # Check if weights are provided and normalize them
  if (!is.null(weights)) {
    # Ensure weights are normalized
    normalized_weights <- normalize_weights(weights)
    # Calculate weighted sum of z-scores
    weighted_sum <- weighted_sum_z_scores(z_scores, normalized_weights)
    # Calculate combined z-score with weights
    combined_z_score <- weighted_sum / sqrt(sum(normalized_weights ^ 2))
  } else {
    # Calculate combined z-score without weights
    combined_z_score <- sum(z_scores) / sqrt(length(z_scores))
  }

  return(combined_z_score)
}

#' Normalize Weights
#'
#' Normalizes a given vector of weights.
#'
#' @param weights A numeric vector of weights.
#'
#' @return A numeric vector of normalized weights.
normalize_weights <- function(weights) {
  weights / sum(weights)
}

#' Calculate Weighted Sum of Z-Scores
#'
#' Calculates the weighted sum of z-scores given z-scores and their
#' corresponding weights.
#'
#' @param z_scores A numeric vector of z-scores.
#' @param weights  A numeric vector of weights for each z-score.
#'
#' @return         Weighted sum of the z-scores.
weighted_sum_z_scores <- function(z_scores, weights) {
  sum(z_scores * weights)
}

### $$$$$$$$$$$$$ Some new $$$$$$$$$$$$$ ###

# ========================================================
# Load and preprocess functions
# ========================================================

#' Load Patient Data into Seurat Object
#'
#' This function loads data for a given patient from a specified directory,
#' constructs the data path, reads the data using Read10X, and creates a Seurat
#' object with metadata.
#'
#' @param patient              List containing patient ID and type.
#' @param base_path            Base path to the data directory.
#' @param analysis_prefix      Prefix for the analysis directory.
#' @param count_default_suffix Suffix for the count directory.
#' @param output_folder_suffix Suffix for the output folder.
#' @param feature_matrix_dir   Directory name for the feature matrix.
#'
#' @return                     A Seurat object with loaded data and metadata.
load_into_seurat <- function(patient, base_path, analysis_prefix,
                             count_default_suffix, output_folder_suffix,
                             feature_matrix_dir) {
  patient_id <- patient$id
  patient_type <- patient$type

  data_dir <- construct_data_dir(base_path, patient_id, analysis_prefix,
                                 count_default_suffix, output_folder_suffix,
                                 feature_matrix_dir)
  data <- Read10X(data.dir = data_dir)

  seurat_object <- create_seurat_object(data, patient_id, patient_type)
  return(seurat_object)
}

#' Construct Data Directory Path
#'
#' This function constructs the data directory path for a given patient based
#' on specified parameters.
#'
#' @param base_path            Base path to the data directory.
#' @param patient_id           Patient ID.
#' @param analysis_prefix      Prefix for the analysis directory.
#' @param count_default_suffix Suffix for the count directory.
#' @param output_folder_suffix Suffix for the output folder.
#' @param feature_matrix_dir   Directory name for the feature matrix.
#'
#' @return                     The constructed data directory path.
construct_data_dir <- function(base_path, patient_id, analysis_prefix,
                               count_default_suffix, output_folder_suffix,
                               feature_matrix_dir) {
  file.path(base_path, patient_id, "analysis",
            paste0(analysis_prefix, patient_id, count_default_suffix),
            paste0(patient_id, output_folder_suffix), feature_matrix_dir)
}

#' Create Seurat Object and Add Metadata
#'
#' This function creates a Seurat object from the given data and adds metadata
#' including patient ID and type.
#'
#' @param data         Data to be loaded into the Seurat object.
#' @param patient_id   Patient ID.
#' @param patient_type Patient type (e.g., Early, Late).
#'
#' @return             A Seurat object with loaded data and metadata.
create_seurat_object <- function(data, patient_id, patient_type) {
  seurat_object <- CreateSeuratObject(counts = data, min.features = 200,
                                      min.cells = 50)
  seurat_object <- RenameCells(seurat_object, add.cell.id = patient_id)
  seurat_object$patient <- patient_id
  seurat_object$type <- patient_type
  return(seurat_object)
}

#' Preprocess Seurat Object
#'
#' This function preprocesses a Seurat object by calculating the percentage of
#' mitochondrial genes, filtering cells, normalizing data using SCTransform,
#' and annotating cells using SingleR.
#'
#' @param seurat_object   A Seurat object to be preprocessed.
#' @param verbose         Boolean indicating whether to print detailed messages.
#' @param blueprint_encode Reference data for SingleR annotation.
#'
#' @return               A preprocessed Seurat object.
preprocess_seurat <- function(seurat_object, verbose = TRUE, blueprint_encode) {
  seurat_object <- calculate_percent_mt(seurat_object)
  seurat_object <- filter_cells(seurat_object)
  seurat_object <- normalize_data(seurat_object, verbose = verbose)
  seurat_object <- annotate_cells_with_singler(seurat_object, blueprint_encode)

  return(seurat_object)
}

#' Calculate Percentage of Mitochondrial Genes
#'
#' This function calculates the percentage of mitochondrial genes for each cell
#' in a Seurat object.
#'
#' @param seurat_object A Seurat object.
#'
#' @return              A Seurat object with added mitochondrial gene
#'                      percentage metadata.
calculate_percent_mt <- function(seurat_object) {
  mitochondrial_genes <- grep("^MT-", rownames(seurat_object), value = TRUE)
  seurat_object[["percent.mt"]] <-
    PercentageFeatureSet(seurat_object, features = mitochondrial_genes)
  return(seurat_object)
}

#' Filter Cells Based on Mitochondrial Content and RNA Count
#'
#' This function filters cells in a Seurat object based on mitochondrial
#' content and RNA count thresholds.
#'
#' @param seurat_object A Seurat object.
#' @param mt_threshold  Threshold for mitochondrial gene percentage.
#' @param min_rna       Minimum RNA count for cells to be retained.
#' @param max_rna       Maximum RNA count for cells to be retained.
#'
#' @return              A filtered Seurat object.
filter_cells <- function(seurat_object, mt_threshold = 25, min_rna = 1000,
                         max_rna = 15000) {
  seurat_object <-
    subset(seurat_object,
           subset = percent.mt < mt_threshold &
             nCount_RNA > min_rna & nCount_RNA < max_rna)
  return(seurat_object)
}

#' Normalize and Stabilize Variance Using SCTransform
#'
#' This function normalizes and stabilizes variance in a Seurat object using
#' the SCTransform method.
#'
#' @param seurat_object A Seurat object.
#' @param verbose       Boolean indicating whether to print detailed messages.
#'
#' @return              A normalized Seurat object.
#' @note                Consider adding `residual_type = "pearson"` to
#'                      SCTransform if corrected UMI counts are needed.
normalize_data <- function(seurat_object, verbose = FALSE) {
  seurat_object <-
    SCTransform(seurat_object, vars.to.regress = c("nCount_RNA", "percent.mt"),
                return.only.var.genes = FALSE, verbose = verbose,
                conserve.memory = TRUE)
  return(seurat_object)
}

#' Annotate Cells Using SingleR
#'
#' This function annotates cells in a Seurat object using SingleR and adds
#' blueprint labels and p-values to the Seurat object.
#'
#' @param seurat_object    A Seurat object.
#' @param blueprint_encode Reference data for SingleR annotation.
#'
#' @return                 A Seurat object with added blueprint labels and
#'                         p-values.
annotate_cells_with_singler <- function(seurat_object, blueprint_encode) {
  singler_results <- SingleR(
    test = seurat_object[["SCT"]]@counts,
    ref = blueprint_encode,
    labels = blueprint_encode$label.main
  )

  seurat_object$blueprint_labels <- singler_results$labels
  seurat_object$blueprint_pvals <- singler_results$scores

  return(seurat_object)
}

# ========================================================
# Data Integration and Batch Correction functions
# ========================================================

# Function to check if Seurat objects are ready for integration
is_seurat_ready_integration <- function(seurat_list, patient_list) {
  invisible(lapply(seq_along(seurat_list), function(i) {
    seurat_object <- seurat_list[[i]]
    patient_id <- patient_list[[i]]$id
    if ("RNA" %in% names(seurat_object@assays) &&
          ncol(GetAssayData(seurat_object,
                            assay = "RNA",
                            layer = "counts")) > 0) {

      message(paste("Seurat object for patient", patient_id,
                    "is ready for integration."))
    } else {
      stop(paste("Seurat object for patient", patient_id,
                 "is not ready for integration."))
    }
  }))
}

# ========================================================
# Clustering functions
# ========================================================

# Function to find the best clustering resolution using silhouette scores
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

# Function to plot silhouette scores
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

# Helper function to set clusters based on the best resolution
set_best_clusters <- function(seurat_obj, best_resolution) {
  # Assign clusters based on the best resolution
  cluster_column_name <- paste("integrated_snn_res.", best_resolution, sep = "")
  seurat_obj$seurat_clusters <- seurat_obj@meta.data[, cluster_column_name]

  # Set the active identity class to the new clusters
  Idents(seurat_obj) <- "seurat_clusters"

  return(seurat_obj)
}

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
    group_by(cluster) %>%
    top_n(n = n_top_genes, wt = avg_log2FC)

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

# ========================================================
# Metacell matrix generation functions
# ========================================================

#' Generate Metacell Matrices for ARACNe Analysis
#'
#' @param seurat_obj   A Seurat object with integrated data.
#' @param out_dir      Output directory to save the metacell matrices.
#' @param out_name     Prefix for the saved metacell matrix files.
#' @param size_thresh  Minimum cluster size to consider.
#' @param num_neighbors Number of neighbors to use for each metacell.
#' @param sub_size     Target number of cells in each metacell.
#'
#' @return A list of metacell matrices, one per cluster.
generate_metacell_matrices <-
  function(seurat_obj, out_dir, out_name, size_thresh = 50, num_neighbors = 5,
           sub_size = 200) {
    counts_matrix <- as.matrix(seurat_obj[["SCT"]]@counts)
    clustering <- seurat_obj$seurat_clusters

    meta_mats <-
      create_metacell_matrices(counts_matrix, clustering, num_neighbors,
                               sub_size, out_dir, out_name, size_thresh)

    return(meta_mats)
  }

#' Create Metacell Matrices for ARACNe Analysis
#'
#' This function takes a gene expression matrix and its corresponding
#' clustering, creates metacell matrices for each cluster, and saves these
#' matrices. It is designed for use with ARACNe to facilitate analysis of gene
#' regulatory networks.
#'
#' @param dat_mat       Matrix of raw gene expression (genes X samples).
#' @param num_neighbors Number of neighbors to use for each metacell.
#' @param clustering    Vector of cluster labels for each sample in `dat_mat`.
#' @param sub_size      Target number of cells in each metacell for ARACNe
#'                      analysis.
#' @param out_dir       Directory where metacell matrices will be saved.
#' @param out_name      Prefix for saved metacell matrix files.
#' @param size_thresh   Minimum cluster size; clusters smaller than this will
#'                      be ignored.
#'
#' @return              A list of metacell matrices, one per cluster.
create_metacell_matrices <-
  function(dat_mat, clustering, num_neighbors = 10, sub_size = 200,
           out_dir, out_name = "", size_thresh = 50) {
    if (is.null(dat_mat) || ncol(dat_mat) == 0) {
      stop("Input data matrix is empty or NULL.")
    }

    if (length(unique(clustering)) <= 1) {
      stop("Insufficient unique clusters for processing.")
    }

    # Generate cluster-specific matrices and filter out empty ones
    clust_mats <- generate_cluster_matrices(dat_mat, clustering, size_thresh)
    clust_mats <- filter_non_empty_matrices(clust_mats)

    # Initialize list to store metacell matrices
    meta_mats <- list()

    if (length(clust_mats) > 0) {
      for (i in seq_along(clust_mats)) {
        meta_mat <- process_cluster(clust_mats[[i]], num_neighbors, i, out_dir,
                                    out_name, sub_size)
        meta_mats[[i]] <- meta_mat
      }
    } else {
      cat("No valid cluster matrices to process.\n")
    }

    return(meta_mats)
  }

#' Generate Cluster-Specific Matrices
#'
#' This function divides a gene expression matrix into sub-matrices based on
#' cluster labels, ensuring each sub-matrix contains only the data for a
#' specific cluster. Clusters with a size below the specified threshold are
#' ignored.
#'
#' @param dat_mat     Matrix of raw gene expression data, with genes as rows
#'                    and samples as columns.
#' @param clustering  A vector of cluster labels corresponding to each column
#'                    in `dat_mat`.
#' @param size_thresh Minimum size of clusters to be considered. Clusters
#'                    smaller than this threshold will be ignored.
#'
#' @return            A list of matrices, each representing gene expression
#'                    data for a specific cluster.
generate_cluster_matrices <- function(dat_mat, clustering, size_thresh) {
  cluster_matrices(dat_mat, clustering, size_thresh = size_thresh)
}

#' Generate and Optionally Save Cluster-Specific Matrices
#'
#' Splits the data matrix into cluster-specific matrices based on provided
#' cluster labels. Can save the resulting matrices to files if a save path is
#' specified.
#'
#' @param dat_mat     Data matrix to be split (features x samples).
#' @param clust       Clustering labels for samples.
#' @param save_path   Optional path for saving the resulting matrices; if not
#'                    provided, matrices are returned in a list.
#' @param save_pref   Optional prefix for file names when saving matrices.
#' @param size_thresh Minimum number of samples required for a cluster to be
#'                    processed; default is 100.
#'
#' @return            A list of matrices, one for each cluster, if `save_path`
#'                    is not provided. Otherwise, nothing is explicitly
#'                    returned.
cluster_matrices <- function(dat_mat, clust, save_path = NA,
                             save_pref = "", size_thresh = 100) {
  clust_table <- table(clust)
  clust_mats <- list()

  for (i in seq_along(clust_table)) {
    if (clust_table[i] > size_thresh) {
      clust_cells <- which(clust == names(clust_table)[i])

      # Check if clust_cells contains valid indices
      if (length(clust_cells) == 0 || max(clust_cells) > ncol(dat_mat)) {
        next  # Skip to the next iteration if no valid cells or out of bounds
      }

      clust_mat <- dat_mat[, clust_cells, drop = FALSE]
      clust_mat <- clust_mat[rowSums(clust_mat) >= 1, , drop = FALSE]

      if (is.na(save_path)) {
        clust_mats[[names(clust_table)[i]]] <- clust_mat
      } else {
        save_cluster_matrix(clust_mat, save_path,
                            save_pref, names(clust_table)[i])
      }
    }
  }

  if (is.na(save_path)) {
    return(clust_mats)
  }
}

#' Save a Cluster Matrix to an RDS File
#'
#' Saves a given cluster matrix to an RDS file, constructing the file name
#' based on provided parameters.
#'
#' @param clust_mat  The cluster-specific matrix to save.
#' @param save_path  Directory path where the file will be saved.
#' @param save_pref  Prefix to be added to the file name.
#' @param clust_name Name of the cluster, used in the file name.
save_cluster_matrix <- function(clust_mat, save_path, save_pref, clust_name) {
  file_path <- file.path(save_path, paste0(save_pref, "_", clust_name, ".rds"))
  saveRDS(clust_mat, file = file_path)
}

#' Filter Out Empty Cluster Matrices
#'
#' Removes any null entries from a list of matrices. This is typically used to
#' exclude cluster-specific matrices that might have been deemed too small or
#' otherwise invalid.
#'
#' @param clust_mats A list of matrices, where each matrix corresponds to a
#'                   cluster's gene expression data.
#'
#' @return           A filtered list of matrices, with null entries removed.
filter_non_empty_matrices <- function(clust_mats) {
  clust_mats <- Filter(function(x) !is.null(x) && nrow(x) > 0, clust_mats)
  return(clust_mats)
}

#' Process Each Cluster to Generate Metacell Matrix
#'
#' For a given cluster's gene expression matrix, this function generates a
#' metacell matrix by considering the specified number of neighbors. It saves
#' two versions of the metacell matrix: one with all cells and another with a
#' subset (if the original exceeds the `sub_size`). Both matrices are saved to
#' files.
#'
#' @param mat           A matrix representing the gene expression data for a
#'                      single cluster.
#' @param num_neighbors The number of neighbors to consider for each metacell.
#' @param cluster_idx   The index of the current cluster being processed.
#' @param out_dir       The directory where output files will be saved.
#' @param out_name      A prefix to be added to the names of the output files.
#' @param sub_size      The maximum number of cells to include in the subsetted
#'                      metacell matrix.
#'
#' @return              A metacell matrix for the cluster, potentially
#'                      subsetted and transformed.
process_cluster <- function(mat, num_neighbors, cluster_idx, out_dir,
                            out_name, sub_size) {
  # Generate metacell matrix
  meta_mat <- meta_cells(mat, num_neighbors)

  # Save the complete metacell matrix
  save_meta_mat(meta_mat, out_dir,
                paste0(out_name, "_clust-", cluster_idx, "-metaCells_all"),
                subset = FALSE)

  # Subset if necessary and apply CPM transformation
  if (sub_size < ncol(meta_mat)) {
    meta_mat <- meta_mat[, sample(colnames(meta_mat), sub_size)]
  }
  meta_mat <- cpmt_transform(meta_mat)

  # Save the subsetted and transformed metacell matrix
  save_meta_mat(meta_mat, out_dir,
                paste0(out_name, "_clust-", cluster_idx, "-metaCells"),
                subset = TRUE)

  return(meta_mat)
}

#' Perform Counts Per Million (CPM) Normalization
#'
#' This function normalizes gene expression data using Counts Per Million (CPM)
#' normalization, with an option for subsequent log2 transformation.
#'
#' @param dat_mat Matrix of gene expression data, with genes as rows and samples
#'                as columns.
#' @param l2      Logical indicating whether to apply log2 transformation after
#'                CPM normalization.
#'
#' @return        A matrix with CPM-normalized (and optionally
#'                log2-transformed) values.
cpmt_transform <- function(dat_mat, l2 = FALSE) {
  # Calculate CPM
  cpm_mat <- calculate_cpm(dat_mat)

  # Apply log2 transformation if specified
  if (l2) {
    cpm_mat <- log2_transform(cpm_mat)
  }

  return(cpm_mat)
}

#' Calculate Counts Per Million (CPM)
#'
#' Converts raw counts to CPM for normalization across samples.
#'
#' @param dat_mat Raw counts matrix with genes as rows and samples as columns.
#'
#' @return        CPM-normalized matrix.
calculate_cpm <- function(dat_mat) {
  t(t(dat_mat) / colSums(dat_mat) * 1e6)
}

#' Apply Log2 Transformation
#'
#' Transforms CPM-normalized values using log2, adding 1 to avoid log of zero.
#'
#' @param cpm_mat Matrix of CPM-normalized gene expression values.
#'
#' @return        Log2-transformed matrix.
log2_transform <- function(cpm_mat) {
  log2(cpm_mat + 1)
}

#' Save Metacell Matrix to File
#'
#' Saves a given metacell matrix to a file, constructing the file name from the
#' provided directory, file prefix, and an indicator of whether the matrix has
#' been subsetted.
#'
#' @param meta_mat    The metacell matrix to be saved.
#' @param out_dir     The directory where the file will be saved.
#' @param file_prefix The prefix to be used in constructing the file name.
#' @param subset      A boolean flag indicating whether the matrix is a
#'                    subsetted version.
#'
#' @return            None; the function's primary effect is to write a file
#'                    to disk.
save_meta_mat <- function(meta_mat, out_dir, file_prefix, subset) {
  # Ensure the output directory exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Construct the file name with appropriate suffix based on the subset flag
  file_suffix <- ifelse(subset, "_sub", "_all")
  file_name <- paste0(out_dir, "/", file_prefix, file_suffix, ".txt")

  # Call aracne_table to save the data appropriately
  aracne_table(meta_mat, file_name, subset)
}

#' Save Data Matrix for ARACNe Analysis
#'
#' Formats and saves a gene expression matrix for use as input to ARACNe,
#' optionally subsetting the matrix to a maximum number of samples for
#' efficiency.
#'
#' @param dat_mat  A matrix of data with genes as rows and samples as columns.
#' @param out_file Path and base name for the output file(s).
#' @param subset   Logical indicating whether to subset the matrix to 500
#'                 samples.
aracne_table <- function(dat_mat, out_file, subset = TRUE) {
  # Remove duplicate genes
  dat_mat <- remove_duplicate_genes(dat_mat)

  # Save the full or subsetted matrix as an RDS file
  save_matrix_rds(dat_mat, out_file, subset)

  # Prepare and save the matrix in TSV format
  save_matrix_tsv(dat_mat, out_file, subset)
}

#' Remove Duplicate Genes from the Matrix
#'
#' @param dat_mat A matrix with genes as rows and samples as columns.
#'
#' @return        Matrix with duplicate genes removed.
remove_duplicate_genes <- function(dat_mat) {
  dat_mat[!duplicated(rownames(dat_mat)), ]
}

#' Save Matrix as RDS File
#'
#' @param dat_mat  A matrix with genes as rows and samples as columns.
#' @param out_file Base path and name for the output RDS file.
#' @param subset   Logical indicating if the matrix should be subsetted.
save_matrix_rds <- function(dat_mat, out_file, subset) {
  if (subset) {
    dat_mat <- subset_matrix_samples(dat_mat, 500)
  }
  saveRDS(dat_mat, file = paste0(out_file, ".rds"))
}

#' Subset Matrix to a Specific Number of Samples
#'
#' @param dat_mat     A matrix with genes as rows and samples as columns.
#' @param max_samples The maximum number of samples to include in the subset.
#'
#' @return            A subsetted matrix.
subset_matrix_samples <- function(dat_mat, max_samples) {
  dat_mat[, sample(colnames(dat_mat), min(ncol(dat_mat), max_samples))]
}

#' Save Matrix as TSV File
#'
#' @param dat_mat  A matrix with genes as rows and samples as columns.
#' @param out_file Base path and name for the output TSV file.
#' @param subset   Logical indicating if the matrix should be subsetted before
#'                 saving.
save_matrix_tsv <- function(dat_mat, out_file, subset) {
  if (subset) {
    dat_mat <- subset_matrix_samples(dat_mat, 500)
  }
  formatted_matrix <- format_matrix_for_tsv(dat_mat)
  write.table(formatted_matrix, file = paste0(out_file, ".tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

#' Format Matrix for Saving as TSV
#'
#' @param dat_mat A matrix with genes as rows and samples as columns.
#'
#' @return        A matrix formatted for saving as TSV, including header row.
format_matrix_for_tsv <- function(dat_mat) {
  sample_names <- colnames(dat_mat)
  gene_ids <- rownames(dat_mat)
  rbind(c("gene", sample_names), cbind(gene_ids, dat_mat))
}

#' Generate Meta Cell Matrix
#'
#' Creates a meta cell matrix by aggregating information from each cell's
#' nearest neighbors. This can be useful for imputing missing values or
#' enhancing signal in sparse datasets.
#'
#' @param dat_mat       A matrix of raw gene expression data (genes x samples).
#' @param num_neighbors The number of nearest neighbors to consider for each
#'                      cell.
#' @param sub_size      Optional; if specified, subsets the resulting meta cell
#'                      matrix to this number of cells.
#'
#' @return              A matrix representing meta cells, potentially
#'                      subsetted.
meta_cells <- function(dat_mat, num_neighbors = 10, sub_size = NA) {
  if (num_neighbors >= ncol(dat_mat)) {
    stop("num_neighbors must be less than the number of columns in dat_mat")
  }
  # Compute distance matrix based on Pearson correlation
  dist_mat <- compute_distance_matrix(dat_mat)

  # Identify nearest neighbors for each sample
  knn_neighbors <- find_knn(dist_mat, num_neighbors)

  # Create imputed matrix based on nearest neighbors
  imp_mat <- impute_matrix(dat_mat, knn_neighbors)

  # Subset the imputed matrix if sub_size is specified and valid
  if (!is.na(sub_size) && sub_size > 0 && sub_size <= ncol(imp_mat)) {
    imp_mat <- subset_matrix(imp_mat, sub_size)
  }

  return(imp_mat)
}

#' Find K-Nearest Neighbors
#'
#' Identifies the k-nearest neighbors for each sample based on the distance
#' matrix.
#'
#' @param dist_mat A distance matrix.
#' @param k        The number of neighbors to identify.
#'
#' @return         A matrix indicating the indices of k-nearest neighbors for
#'                 each sample.
find_knn <- function(dist_mat, k) {
  # Ensure the distance matrix is in matrix format
  dist_mat <- as.matrix(dist_mat)

  # Apply over each row to find the indices of the k-nearest neighbors
  knn_indices <- t(apply(dist_mat, 1, function(x) {
    actual_k <- min(k, length(x) - 1)
    order(x)[2:(actual_k + 1)]
  }))

  return(knn_indices)
}

#' Impute Matrix
#'
#' Creates an imputed matrix by aggregating the expression of each sample with
#' its k-nearest neighbors.
#'
#' @param dat_mat       A matrix of gene expression data (genes x samples).
#' @param knn_neighbors A matrix of k-nearest neighbor indices for each sample.
#'
#' @return              An imputed gene expression matrix.
impute_matrix <- function(dat_mat, knn_neighbors) {
  # Initialize the imputed matrix
  imp_mat <- matrix(0, nrow = nrow(dat_mat), ncol = ncol(dat_mat))
  colnames(imp_mat) <- colnames(dat_mat)
  rownames(imp_mat) <- rownames(dat_mat)

  # Iterate over each sample to impute based on nearest neighbors
  for (i in seq_len(ncol(dat_mat))) {
    # Retrieve neighbor indices for the current sample
    neighbor_cols <- c(i, knn_neighbors[i, ])

    # Ensure all referenced indices are within bounds
    if (any(neighbor_cols > ncol(dat_mat))) {
      stop(paste("Out of bounds error at sample", i,
                 ": Neighbor indices",
                 toString(neighbor_cols[neighbor_cols > ncol(dat_mat)]),
                 "are greater than the number of columns", ncol(dat_mat)))
    }

    # Aggregate data from the original matrix using the neighbor indices
    imp_mat[, i] <-
      rowSums(dat_mat[, neighbor_cols, drop = FALSE], na.rm = TRUE)
  }

  return(imp_mat)
}

#' Subset Matrix
#'
#' Subsets a matrix to a specified number of columns (samples), chosen
#' randomly.
#'
#' @param mat      A matrix to be subsetted.
#' @param sub_size The number of columns to retain in the subset.
#'
#' @return         A subsetted matrix.
subset_matrix <- function(mat, sub_size) {
  mat[, sample(ncol(mat), sub_size)]
}

#' Plot Cluster Frequencies by Treatment
#'
#' @param seurat_obj       A Seurat object with integrated data.
#' @param plot_output_path Directory to save the plot.
#' @param col_names        Vector of column names for the combined data frame.
#' @param plot_title       Title of the plot.
#' @param binwidth         Width of bins in the dot plot.
#'
#' @return None. The function saves a plot to the specified directory.
plot_cluster_freq_by_treatment <-
  function(seurat_obj, plot_output_path,
           col_names = c("Early_p1", "Early_p2", "Early_p3", "Late_p1",
                         "Late_p2", "Late_p3"),
           plot_title = "Cluster Frequency by Early vs Late",
           binwidth = 0.01) {

    cluster_freq_table <-
      table(seurat_obj$patient, seurat_obj$seurat_clusters, seurat_obj$type)

    early_data <- calculate_cluster_frequencies(cluster_freq_table, 1)
    late_data <- calculate_cluster_frequencies(cluster_freq_table, 2)

    combined_data <-
      combine_cluster_frequencies(early_data, late_data, col_names)

    plot_data <- melt_cluster_frequencies(combined_data)

    plot_cluster_frequencies(plot_data, plot_title,
                             file.path(plot_output_path,
                                       "cluster_frequencies_by_treatment.png"),
                             binwidth)
  }

#' Calculate Cluster Frequencies
#'
#' @param cluster_freq_table Table of cluster frequencies.
#' @param treatment_index Index of the treatment type (1 for early, 2 for late).
#'
#' @return A data frame of normalized cluster frequencies.
calculate_cluster_frequencies <- function(cluster_freq_table, treatment_index) {
  freq_data <- as.data.frame.matrix(cluster_freq_table[, , treatment_index])
  freq_data <- freq_data[rowSums(freq_data) > 0, ]
  freq_data <- apply(freq_data, 1, function(x) {
    x / sum(x)
  })
  return(freq_data)
}

#' Combine Cluster Frequencies
#'
#' @param early_data Data frame of early treatment cluster frequencies.
#' @param late_data Data frame of late treatment cluster frequencies.
#' @param col_names Vector of column names for the combined data frame.
#'
#' @return A combined data frame of early and late treatment cluster
#'         frequencies.
combine_cluster_frequencies <- function(early_data, late_data, col_names) {
  combined <- cbind(early_data, late_data)
  colnames(combined) <- col_names
  return(combined)
}

#' Melt Cluster Frequencies
#'
#' @param combined_data Combined data frame of cluster frequencies.
#'
#' @return A melted data frame suitable for plotting.
melt_cluster_frequencies <- function(combined_data) {
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
}

#' Plot Cluster Frequencies
#'
#' @param plot_data Data frame of cluster frequencies to be plotted.
#' @param plot_title Title of the plot.
#' @param output_path File path to save the plot.
#' @param binwidth Width of bins in the dot plot.
#'
#' @return None. The function saves a plot to the specified file path.
plot_cluster_frequencies <-
  function(plot_data, plot_title, output_path, binwidth) {
    theme_update(plot.title = element_text(hjust = 0.5))

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
  }

# ========================================================
# Running ARACNe Analysis functions
# ========================================================

#' Prepare and Save Expression Matrix for ARACNe Analysis
#'
#' This function prepares and saves an expression matrix for use in ARACNe
#' analysis. It first removes duplicate genes, then subsets the matrix if
#' necessary, and finally saves the matrix in both RDS and TSV formats.
#'
#' @param base_output_path  The base path for the output files.
#' @param file_suffix       The suffix to append to the input file name.
#'
prep_and_save_expr_for_aracne <- function(base_output_path, file_suffix) {
  files <- list.files(base_output_path, pattern = paste0(file_suffix, "$"),
                      full.names = TRUE)
  expr_files <- lapply(files, function(file_path) {
    expr_data <- read.table(file_path, header = TRUE, sep = "\t",
                            row.names = 1)
    out_path <- gsub(".tsv$", "_for_aracne.tsv", file_path)
    save_matrix_for_aracne(expr_data, out_path)
    return(out_path)
  })
  return(expr_files)
}

#' Save Expression Matrix for ARACNe
#'
#' Saves an expression matrix in a format required by ARACNe3.
#'
#' @param expression_matrix The normalized expression matrix to be saved.
#' @param output_file Path for the output file.
save_matrix_for_aracne <- function(expression_matrix, output_file) {
  # Ensure the matrix has proper row and column names
  if (is.null(colnames(expression_matrix)) ||
        is.null(rownames(expression_matrix))) {
    stop("Expression matrix must have row and column names.")
  }

  # Ensure there is an extra column in the header (if not already present)
  cat("\t", file = output_file, append = FALSE)
  suppressWarnings(
    write.table(expression_matrix, file = output_file, sep = "\t",
                quote = FALSE, row.names = TRUE, col.names = TRUE,
                append = TRUE)
  )
}

# Function to run ARACNe for each expression and regulator file
run_aracne_for_all <- function(aracne_bin, expression_files, regulator_files,
                               output_base_dir, threads, seed) {
  for (reg_name in names(regulator_files)) {
    regulator_file <- regulator_files[[reg_name]]
    for (exp_file in expression_files) {
      exp_file_base <-
        gsub("_all_all.txt_for_aracne.tsv", "", basename(exp_file))
      output_dir <-
        file.path(output_base_dir, paste0(reg_name, "_", exp_file_base))
      create_directories(list(output_dir))
      run_aracne(aracne_bin, exp_file, regulator_file, output_dir, threads,
                 seed)
    }
  }
}

#' Run ARACNe3
#'
#' Executes ARACNe3 on the provided expression matrix and regulator list.
#'
#' @param aracne_bin Path to the ARACNe3 binary.
#' @param exp_file Path to the expression matrix file.
#' @param regulators_file Path to the file containing regulator gene names.
#' @param output_dir Directory to store ARACNe output.
#' @param threads Number of threads to use for ARACNe computation.
#' @param seed Seed for random number generation in ARACNe.
run_aracne <- function(aracne_bin, exp_file, regulators_file, output_dir,
                       threads = 1, seed = 123) {
  cmd <-
    sprintf("%s -e %s -r %s -o %s --threads %d --seed %d",
            aracne_bin, exp_file, regulators_file, output_dir, threads, seed)
  system(cmd)
}

# ========================================================
# Running VIPER Analysis functions
# ========================================================

# Function to process ARACNe output files
generate_regulon_objects <- function(aracne_output_base_dir, exp_mat,
                                     regulator_files, output_base_path) {
  aracne_output_files <-
    get_all_aracne_files(aracne_output_base_dir, regulator_files)

  if (length(aracne_output_files) == 0) {
    stop("No ARACNe output files found in the directory.")
  }

  regulon_list <- lapply(aracne_output_files, function(aracne_file) {
    aracne_data_for_viper <- prep_aracne_output_for_viper(aracne_file)
    regulon <-
      generate_regulon(aracne_data_for_viper, exp_mat, output_base_path,
                       paste0(basename(aracne_file), "_"))
    return(regulon)
  })

  return(regulon_list)
}

# Helper function to aggregate ARACNe output files from all directories
get_all_aracne_files <- function(base_dir, regulator_files) {
  all_files <- unlist(lapply(names(regulator_files), function(reg_name) {
    reg_dir <- file.path(base_dir, reg_name)
    list.files(reg_dir,
               pattern = "consolidated-net_.*\\.tsv$", full.names = TRUE)
  }))
  return(all_files)
}

# Helper function to load and process ARACNe output file
prep_aracne_output_for_viper <- function(aracne_file) {
  cat("Loading ARACNe output file:", aracne_file, "\n")

  # Load ARACNe output file without headers
  aracne_data <- read.table(aracne_file, header = FALSE, sep = "\t",
                            check.names = FALSE, stringsAsFactors = FALSE,
                            skip = 1)

  # Only include the first three columns
  aracne_data <- aracne_data[, 1:3]

  # Define column names manually
  colnames(aracne_data) <- c("regulator", "target", "mi")

  # Convert the 'mi' column to numeric
  aracne_data$mi <- as.numeric(aracne_data$mi)

  return(aracne_data)
}

# Helper function to generate regulon object from ARACNe output
generate_regulon <- function(aracne_data, exp_mat, output_base_path,
                             file_prefix) {
  # Process ARACNe results for VIPER analysis
  reg_process(aracne_data, exp_mat, output_base_path, file_prefix)

  # Load the pruned regulon object
  pruned_regulon_file <-
    file.path(output_base_path, paste0(file_prefix, "pruned.rds"))
  pruned_regulon <- readRDS(pruned_regulon_file)

  return(pruned_regulon)
}

# Function to run VIPER on a list of regulon objects
run_viper <- function(exp_mat, regulon_list) {
  viper_results <- lapply(regulon_list, function(regulon) {
    viper_scores <- run_viper_analysis(exp_mat, regulon)
    return(viper_scores)
  })

  # Check contents of viper_results
  if (length(viper_results) == 0 || any(sapply(viper_results, is.null))) {
    stop("VIPER results are empty or not properly formed.")
  }

  return(viper_results)
}

# Function to run VIPER analysis on a single regulon
run_viper_analysis <- function(exp_mat, regulon) {
  viper_scores <- tryCatch({
    viper(exp_mat, regulon)
  }, error = function(e) {
    cat("Error in VIPER analysis:", e$message, "\n")
    NULL
  })

  return(viper_scores)
}

# Function to save VIPER results
save_viper_results <- function(viper_results, output_path) {
  saveRDS(viper_results, file = output_path)
}

# ========================================================
# Utility functions
# ========================================================

# Function to create multiple directories if they don't exist
create_directories <- function(dir_paths) {
  for (dir_path in dir_paths) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
  }
}