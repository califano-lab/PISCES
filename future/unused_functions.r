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