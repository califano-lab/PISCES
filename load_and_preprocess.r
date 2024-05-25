library(Seurat)
library(SingleR)

#' Load and Preprocess Patient Data
#'
#' This function loads data for each patient from specified directories,
#' constructs the data path, reads the data using Read10X, creates a Seurat
#' object with metadata, and preprocesses the Seurat object.
#'
#' @param patients       List containing metadata for each patient
#' @param base_data_path Base path to the data directory.
#' @param patient_data_path Relative path to the data directory within each
#'                          patient directory.
#' @param mt_threshold   Threshold for mitochondrial gene percentage.
#' @param min_rna        Minimum number of RNA features.
#' @param max_rna        Maximum number of RNA features.
#' @param my_verbose     Set to TRUE to display verbose messages during the
#'                       analysis.
#'
#' @return               A list of preprocessed Seurat objects for each patient.
load_and_preprocess <- function(patients, base_data_path, patient_data_path,
                                mt_threshold = 25, min_rna = 1000,
                                max_rna = 15000, my_verbose = FALSE) {
  # Load data for each patient into separate Seurat objects
  patient_seurat_list <- lapply(patients, function(patient) {
    tryCatch({
      message("Loading patient: ", patient$id)
      seurat_obj <- load_into_seurat(patient, base_data_path, patient_data_path)
      message("Completed loading for patient: ", patient$id)
      return(seurat_obj)
    }, error = function(e) {
      stop("Error loading patient: ", patient$id, ": ", e$message)
    })
  })

  message("Finished loading patient data into Seurat objects.")

  blueprint_encode <- BlueprintEncodeData()

  # Preprocess each Seurat object
  patient_seurat_list <- lapply(patient_seurat_list, function(p) {
    patient_id <- unique(p$id)
    tryCatch({
      message("Preprocessing Seurat object for patient: ", patient_id)
      p <- preprocess_seurat(p, blueprint_encode,
                             mt_threshold = mt_threshold, min_rna = min_rna,
                             max_rna = max_rna, my_verbose)
      message("Completed preprocessing Seurat object for patient: ", patient_id)
      return(p)
    }, error = function(e) {
      stop("Error preprocessing Seurat object for patient: ",
           patient_id, ": ", e$message)
    })
  })

  message("Finished preprocessing seurat object(s).")

  return(patient_seurat_list)
}


#' Load Patient Data into Seurat Object
#'
#' This function loads data for a given patient from a specified directory,
#' constructs the data path, reads the data using Read10X, and creates a Seurat
#' object with metadata.
#'
#' @param patient              List containing patient metadata.
#' @param base_path            Base path to the data directory.
#' @param patient_data_path    Relative path to the data directory within each
#'                             patient directory.
#'
#' @return                     A Seurat object with loaded data and metadata.
load_into_seurat <- function(patient, base_path, patient_data_path) {
  patient_id <- patient$id

  data_dir <- file.path(base_path, patient_id, patient_data_path)
  data <- Read10X(data.dir = data_dir)

  seurat_object <- create_seurat_object(data, patient)
  return(seurat_object)
}

#' Create Seurat Object and Add Metadata
#'
#' This function creates a Seurat object from the given data and adds metadata
#' for each field in the patient list.
#'
#' @param data         Data to be loaded into the Seurat object.
#' @param patient      List containing patient metadata.
#'
#' @return             A Seurat object with loaded data and metadata.
create_seurat_object <- function(data, patient) {
  seurat_object <-
    CreateSeuratObject(counts = data, min.features = 200, min.cells = 50)
  seurat_object <- RenameCells(seurat_object, add.cell.id = patient$id)

  # Add each metadata field to the Seurat object
  for (field in names(patient)) {
    seurat_object[[field]] <- patient[[field]]
  }

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
preprocess_seurat <- function(seurat_object, blueprint_encode,
                              mt_threshold = 25, min_rna = 1000,
                              max_rna = 15000, verbose = FALSE) {
  seurat_object <- calculate_percent_mt(seurat_object)
  seurat_object <- filter_cells(seurat_object, mt_threshold, min_rna, max_rna)
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
filter_cells <- function(seurat_object, mt_threshold, min_rna, max_rna) {
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