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
Loader <- R6Class( # nolint
  "Loader",
  public = list(
    patients = NULL,
    base_data_path = NULL,
    patient_data_path = NULL,

    #' Initialize the Loader
    #'
    #' @param patients          A list of patients with their metadata.
    #' @param base_data_path    The base path to the data directory.
    #' @param patient_data_path The relative path to the data directory within
    #'                          each patient directory.
    initialize = function(patients, base_data_path, patient_data_path) {
      self$patients <- patients
      self$base_data_path <- base_data_path
      self$patient_data_path <- patient_data_path
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
        CreateSeuratObject(counts = data, min.features = 200, min.cells = 50)
      seurat_object <- RenameCells(seurat_object, add.cell.id = patient$id)

      # Add each metadata field to the Seurat object
      for (field in names(patient)) {
        seurat_object[[field]] <- patient[[field]]
      }

      return(seurat_object)
    }
  )
)