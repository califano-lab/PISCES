library(Seurat)
library(SingleR)

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
           subset = percent.mt < mt_threshold & # nolint
             nCount_RNA > min_rna & nCount_RNA < max_rna) # nolint
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