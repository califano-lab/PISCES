library(R6)

#' TableGenerator Class for generating different tables
#'
#' @field seurat_obj A Seurat object containing the single-cell data.
#' @field base_output_path A string specifying the base path where output
#'        files will be saved.
TableGenerator <- R6Class( # nolint
  "TableGenerator",
  public = list(
    seurat_obj = NULL,
    base_output_path = NULL,

    #' Initialize the TableGenerator Object
    #'
    #' @param seurat_obj       A Seurat object containing the single-cell data.
    #' @param base_output_path A string specifying the base path where output
    #'                         files will be saved.
    initialize = function(seurat_obj, base_output_path) {
      self$seurat_obj <- seurat_obj
      self$base_output_path <- base_output_path
    },

    #' Generate Cell Count Summary Table
    #'
    #' This method generates a summary table containing the pre-QC and post-QC
    #' cell counts for each patient in the Seurat object. The summary table is
    #' saved as a CSV file in the specified output path.
    #'
    #' @return A data frame containing the summary table with columns for
    #'         Patient, PreQC, and PostQC cell counts.
    generate_cell_count_summary = function() {
      patient_ids <- unique(self$seurat_obj@meta.data$id)
      summary_data <- lapply(patient_ids, function(patient_id) {
        patient_data <- subset(self$seurat_obj, id == patient_id)
        pre_qc_count <- ncol(patient_data)
        post_qc_count <- sum(patient_data@meta.data$nFeature_RNA > 200
                             & patient_data@meta.data$nCount_RNA < 2500)
        data.frame(Patient = patient_id, PreQC = pre_qc_count,
                   PostQC = post_qc_count)
      })
      summary_table <- do.call("rbind", summary_data)
      write.csv(summary_table,
                file = file.path(base_output_path, "cell_counts_summary.csv"),
                row.names = FALSE)

      return(summary_table)
    }
  ),

  #############################################################################
  #                           PRIVATE METHODS                                 #
  #############################################################################
  private = list()
)