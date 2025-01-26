library(R6)
library(Seurat)
library(SeuratWrappers)
library(batchelor)

#' Integrator Class for Integrating Seurat Objects
#'
#' This class provides methods to integrate a list of Seurat objects using
#' either the Seurat or fastMNN method.
#'
#' @field seurat_list A list of Seurat objects to be integrated.
#' @field nfeatures Number of features to select for integration.
#' @field reference Index of the reference dataset for integration.
#' @field integration_type Type of integration method ("seurat" or "fastMNN").
#' @field verbose Boolean indicating whether to print detailed messages.
Integrator <- R6Class( # nolint
  "Integrator",
  public = list(
    seurat_list = NULL,
    nfeatures = 4000,
    reference = 1,
    integration_type = "seurat",
    verbose = FALSE,

    #' Initialize the Integrator Object
    #'
    #' @param seurat_list      A list of Seurat objects to be integrated.
    #' @param nfeatures        Number of features to select for integration.
    #'                         Defaults to 4000.
    #' @param reference        Index of the reference dataset for integration.
    #'                         Defaults to 1.
    #' @param integration_type Type of integration method ("seurat" or
    #'                         "fastMNN"). Defaults to "seurat".
    #' @param verbose          Boolean indicating whether to print detailed
    #'                         messages. Defaults to FALSE.
    initialize = function(seurat_list, nfeatures = 4000, reference = 1,
                          integration_type = "seurat",
                          verbose = FALSE) {
      self$seurat_list <- seurat_list
      self$nfeatures <- nfeatures
      self$reference <- reference
      self$integration_type <- integration_type
      self$verbose <- verbose
    },

    #' Integrate Data
    #'
    #' This function integrates a list of Seurat objects using the specified
    #' integration method ("seurat" or "fastMNN").
    #'
    #' @return An integrated Seurat object.
    integrate_data = function() {
      if (self$integration_type == "seurat") {
        return(private$integrate_data_seurat(self$seurat_list,
                                             self$verbose))
      } else if (self$integration_type == "fastMNN") {
        return(private$integrate_data_fastmnn(self$seurat_list))
      } else {
        stop(paste0("Unsupported integration type.",
                    " Choose either 'seurat' or 'fastMNN'."))
      }
    }
  ),

  #############################################################################
  #                           PRIVATE METHODS                                 #
  #############################################################################
  private = list(
    #' Integrate Data Using Seurat Method
    #'
    #' This function integrates a list of Seurat objects using the Seurat
    #' integration method.
    #'
    #' @param patient_seurat_list A list of Seurat objects to be integrated.
    #' @param verbose             Boolean indicating whether to print detailed
    #'                            messages.
    #'
    #' @return                    An integrated Seurat object.
    integrate_data_seurat = function(patient_seurat_list, verbose) {
      private$is_seurat_ready_integration(patient_seurat_list)

      features_to_integrate <-
        SelectIntegrationFeatures(object.list = patient_seurat_list,
                                  nfeatures = self$nfeatures)

      patient_seurat_list <-
        PrepSCTIntegration(object.list = patient_seurat_list,
                           anchor.features = features_to_integrate,
                           verbose = verbose)

      patient_seurat_list <- lapply(patient_seurat_list, function(seurat_obj) {
        RunPCA(seurat_obj, features = features_to_integrate,
               verbose = verbose)
      })

      anchors <-
        FindIntegrationAnchors(object.list = patient_seurat_list,
                               anchor.features = features_to_integrate,
                               dims = 1:30, normalization.method = "SCT",
                               reduction = "rpca", k.anchor = 20,
                               verbose = verbose,
                               reference = self$reference)

      rm(patient_seurat_list, features_to_integrate)

      integrated_seurat <-
        IntegrateData(anchorset = anchors, normalization.method = "SCT",
                      dims = 1:30, verbose = verbose)

      rm(anchors)

      return(integrated_seurat)
    },

    #' Check if Seurat Objects are Ready for Integration
    #'
    #' It verifies that each Seurat object contains the "RNA" assay with
    #' non-empty counts. If a Seurat object is not ready, the function stops
    #' and returns an error message.
    #'
    #' @param seurat_list A list of Seurat objects to be checked.
    #'
    #' @return            An invisible list of messages indicating which Seurat
    #'                    objects are ready for integration.
    is_seurat_ready_integration = function(seurat_list) {
      invisible(lapply(seurat_list, function(seurat_object) {
        patient_id <- seurat_object@meta.data[["id"]][1]
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
    },

    #' Integrate Data Using fastMNN Method
    #'
    #' This function integrates a list of Seurat objects using the fastMNN
    #' integration method.
    #'
    #' @param patient_seurat_list A list of Seurat objects to be integrated.
    #'
    #' @return                    An integrated Seurat object.
    integrate_data_fastmnn = function(patient_seurat_list) {
      # Normalize and find variable features for each Seurat object
      patient_seurat_list <- lapply(patient_seurat_list, function(seurat_obj) {
        seurat_obj <- NormalizeData(seurat_obj)
        seurat_obj <- FindVariableFeatures(seurat_obj)
        return(seurat_obj)
      })

      # Run FastMNN integration
      integrated_seurat <-
        RunFastMNN(object.list = patient_seurat_list, features = features)

      integrated_seurat <-
        RunUMAP(integrated_seurat, reduction = "mnn", dims = 1:30)
      integrated_seurat <-
        FindNeighbors(integrated_seurat, reduction = "mnn", dims = 1:30)
      integrated_seurat <-
        FindClusters(integrated_seurat)

      return(integrated_seurat)
    }
  )
)
