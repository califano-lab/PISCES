library(R6)

#' Loader Class for Loading Patient Data into Seurat Objects
#'
#' This R6 class provides methods to load patient data from specified
#' directories, create Seurat objects, and add metadata to the Seurat objects.
#'
#' @field patients A list of patients with their metadata.
#' @field base_data_path The base path to the data directory.
#' @field patient_data_path The relative path to the data directory within each
#'        patient directory.
#' @field min_cells Genes must be detected in at least this many cells WITHIN
#'        a patient to be kept.
#' @field min_features Cells must express at least this many genes to be kept.
Loader <- R6Class( # nolint
  "Loader",
  public = list(
    patients = NULL,
    base_data_path = NULL,
    patient_data_path = NULL,
    min_cells = 0,
    min_features = 0,

    #' Initialize the Loader
    #'
    #' @param patients          A list of patients with their metadata.
    #' @param base_data_path    The base path to the data directory.
    #' @param patient_data_path The relative path to the data directory within
    #'                          each patient directory.
    #' @param min_cells         CreateSeuratObject min.cells. Applied PER
    #'                          PATIENT, so it is a much harsher filter for a
    #'                          small patient than a large one. Defaults to 0.
    #' @param min_features      CreateSeuratObject min.features. Defaults to 0.
    #'
    #' @note THESE WERE HARD-CODED AT min.cells = 50, min.features = 200, and
    #'       that combination silently destroys small patients. CreateSeuratObject
    #'       filters GENES FIRST, then cells. Measured on this project's 25
    #'       neutrophil Origins:
    #'
    #'         Origin      cells  genes detected  surviving min.cells = 50
    #'         GSE241184     101            6550                        83
    #'         GSE184198     139            5622                       148
    #'         GSE215403     159            7656                       223
    #'         TS          17217           15517                     10728
    #'
    #'       A 101-cell patient keeps 83 genes, and min.features = 200 then
    #'       requires every cell to express 200 genes that no longer exist - so
    #'       every cell is dropped and the patient becomes an empty object.
    #'
    #'       It also collapses the shared gene space that integration needs:
    #'       across eight of these Origins the intersection is 44 genes at
    #'       min.cells = 50 versus 1,454 at min.cells = 3. SelectIntegrationFeatures
    #'       cannot return more anchors than that intersection allows, so raising
    #'       nfeatures has no effect while this is set high.
    #'
    #'       Defaults are now 0 because gene and cell filtering belong upstream,
    #'       where they can be applied to the whole cohort at once rather than
    #'       per patient. Raise them only for genuinely raw input.
    initialize = function(patients, base_data_path, patient_data_path,
                          min_cells = 0, min_features = 0) {
      self$patients <- patients
      self$base_data_path <- base_data_path
      self$patient_data_path <- patient_data_path
      self$min_cells <- min_cells
      self$min_features <- min_features
    },

    #' Load Data for All Patients
    #'
    #' This method loads data for all patients, converts them into Seurat
    #' objects, and adds metadata.
    #'
    #' @return A list of Seurat objects for each patient.
    load_data = function() {
      patient_seurat_list <- lapply(self$patients, function(patient) {
        tryCatch({
          message("Loading patient: ", patient$id)
          patient_path <- file.path(self$base_data_path, patient$id)
          if (file.exists(paste0(patient_path, ".rds"))) {
            seurat_obj <- private$load_rds_file(patient, patient_path)
          } else {
            seurat_obj <- private$load_bc_matrix(patient)
          }
          message("Completed loading for patient: ", patient$id)
          return(seurat_obj)
        }, error = function(e) {
          stop("Error loading patient: ", patient$id, ": ", e$message)
        })
      })

      message("Finished loading patient data into Seurat objects.")
      return(patient_seurat_list)
    }
  ),

  #############################################################################
  #                           PRIVATE METHODS                                 #
  #############################################################################
  private = list(
    #' Load Patient Data from RDS File into Seurat Object
    #'
    #' This function loads data for a given patient from an RDS file,
    #' reads the data and creates a Seurat object with metadata.
    #'
    #' @param patient              List containing patient metadata.
    #' @param patient_path         Path to the patient .rds file.
    #'
    #' @return                     A Seurat object with loaded data and
    #'                             metadata.
    load_rds_file = function(patient, patient_path) {
      rds_file <- paste0(patient_path, ".rds")
      data <- readRDS(rds_file)
      seurat_object <- private$create_seurat_object(data, patient)
      return(seurat_object)
    },

    #' Load Patient Data from Feature Bc Matrix into Seurat Object
    #'
    #' This function loads data for a given patient from a specified directory,
    #' constructs the data path, reads the data using Read10X, and creates a
    #' Seurat object with metadata.
    #'
    #' @param patient              List containing patient metadata.
    #' @param base_path            Base path to the data directory.
    #' @param patient_data_path    Relative path to the data directory within
    #'                             each patient directory.
    #'
    #' @return                     A Seurat object with loaded data and
    #'                             metadata.
    load_bc_matrix = function(patient) {
      patient_id <- patient$id
      data_dir <-
        file.path(self$base_data_path, patient_id, self$patient_data_path)
      data <- Read10X(data.dir = data_dir)
      seurat_object <- private$create_seurat_object(data, patient)
      return(seurat_object)
    },

    #' Create Seurat Object and Add Metadata
    #'
    #' This function creates a Seurat object from the given data and adds
    #' metadata for each field in the patient list.
    #'
    #' @param data         Data to be loaded into the Seurat object.
    #' @param patient      List containing patient metadata.
    #'
    #' @return             A Seurat object with loaded data and metadata.
    create_seurat_object = function(data, patient) {
      seurat_object <-
        CreateSeuratObject(counts = data,
                           min.features = self$min_features,
                           min.cells = self$min_cells)
      if (ncol(seurat_object) == 0) {
        stop("Patient ", patient$id, " has 0 cells after CreateSeuratObject ",
             "(min.cells = ", self$min_cells, ", min.features = ",
             self$min_features, "). Genes are filtered before cells, so a high ",
             "min.cells on a small patient can leave fewer genes than ",
             "min.features requires.")
      }
      message(sprintf("  %s: %d cells x %d genes", patient$id,
                      ncol(seurat_object), nrow(seurat_object)))
      seurat_object <- RenameCells(seurat_object, add.cell.id = patient$id)

      # Add each metadata field to the Seurat object
      for (field in names(patient)) {
        seurat_object[[field]] <- patient[[field]]
      }

      return(seurat_object)
    }
  )
)