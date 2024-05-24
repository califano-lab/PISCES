library(Seurat)

#' Check if Seurat Objects are Ready for Integration
#'
#' It verifies that each Seurat object contains the "RNA" assay with non-empty
#' counts. If a Seurat object is not ready, the function stops and returns an
#' error message.
#'
#' @param seurat_list A list of Seurat objects to be checked.
#' @param patient_list A list of patient information, where each element
#'                     contains patient details including the 'id'.
#'
#' @return An invisible list of messages indicating which Seurat objects are
#'         ready for integration.
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