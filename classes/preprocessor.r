library(Seurat)
library(SingleR)
library(R6)

#' Preprocessor Class for Preprocessing Seurat Objects
#'
#' This class provides methods to preprocess a list of Seurat objects by
#' calculating mitochondrial gene percentages, filtering cells, normalizing
#' data using SCTransform, and annotating cells using SingleR.
#'
#' @field seurat_list A list of Seurat objects to be preprocessed.
#' @field mt_threshold Threshold for mitochondrial gene percentage.
#' @field min_rna Minimum RNA count for cells to be retained.
#' @field max_rna Maximum RNA count for cells to be retained.
#' @field verbose Boolean indicating whether to print detailed messages.
#' @field blueprint_encode Reference data for SingleR annotation.
Preprocessor <- R6Class( # nolint
  "Preprocessor",
  public = list(
    seurat_list = NULL,
    mt_threshold = 25,
    min_rna = 1000,
    max_rna = 15000,
    verbose = FALSE,
    blueprint_encode = NULL,

    #' Initialize the Preprocessor Object
    #'
    #' @param seurat_list  A list of Seurat objects to be preprocessed.
    #' @param mt_threshold Threshold for mitochondrial gene percentage.
    #'                     Defaults to 25.
    #' @param min_rna      Minimum RNA count for cells to be retained.
    #'                     Defaults to 1000.
    #' @param max_rna      Maximum RNA count for cells to be retained.
    #'                     Defaults to 15000.
    #' @param verbose      Boolean indicating whether to print detailed
    #'                     messages. Defaults to FALSE.
    initialize = function(seurat_list, mt_threshold = 25,
                          min_rna = 1000, max_rna = 15000,
                          verbose = FALSE) {
      self$seurat_list <- seurat_list
      self$mt_threshold <- mt_threshold
      self$min_rna <- min_rna
      self$max_rna <- max_rna
      self$verbose <- verbose
      self$blueprint_encode <- celldex::BlueprintEncodeData()
    },

    #' Preprocess Data
    #'
    #' This function preprocesses a list of Seurat objects by calculating
    #' the percentage of mitochondrial genes, filtering cells, normalizing
    #' data using SCTransform, and annotating cells using SingleR.
    #'
    #' @return A list of preprocessed Seurat objects.
    preprocess_data = function() {
      seurat_list <- lapply(self$seurat_list, function(seurat_object) {
        patient_id <- unique(seurat_object$id)
        tryCatch({
          message("Preprocessing Seurat object for patient: ", patient_id)
          seurat_object <- private$preprocess_seurat(seurat_object)
          message("Completed preprocessing Seurat object for patient: ",
                  patient_id)
          return(seurat_object)
        }, error = function(e) {
          stop("Error preprocessing Seurat object for patient: ",
               patient_id, ": ", e$message)
        })
      })

      message("Finished preprocessing Seurat object(s).")
      return(seurat_list)
    }
  ),

  #############################################################################
  #                           PRIVATE METHODS                                 #
  #############################################################################
  private = list(
    #' Preprocess Seurat Object
    #'
    #' This function preprocesses a Seurat object by calculating the percentage
    #' of mitochondrial genes, filtering cells, normalizing data using
    #' SCTransform, and annotating cells using SingleR.
    #'
    #' @param seurat_object    A Seurat object to be preprocessed.
    #' @param verbose          Boolean indicating whether to print detailed
    #'                         messages.
    #' @param blueprint_encode Reference data for SingleR annotation.
    #'
    #' @return                 A preprocessed Seurat object.
    preprocess_seurat = function(seurat_object) {
      seurat_object <- private$calculate_percent_mt(seurat_object)
      seurat_object <- private$filter_cells(seurat_object)
      seurat_object <- private$normalize_data(seurat_object)
      seurat_object <- private$annotate_cells_with_singler(seurat_object)
      return(seurat_object)
    },

    #' Calculate Percentage of Mitochondrial Genes
    #'
    #' This function calculates the percentage of mitochondrial genes for each
    #' cell in a Seurat object.
    #'
    #' @param seurat_object A Seurat object.
    #'
    #' @return              A Seurat object with added mitochondrial gene
    #'                      percentage metadata.
    calculate_percent_mt = function(seurat_object) {
      mitochondrial_genes <-
        grep("^MT-", rownames(seurat_object), value = TRUE)
      seurat_object[["percent.mt"]] <-
        PercentageFeatureSet(seurat_object, features = mitochondrial_genes)
      return(seurat_object)
    },

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
    filter_cells = function(seurat_object) {
      seurat_object <-
        subset(seurat_object,
               subset = percent.mt < self$mt_threshold &
                 nCount_RNA > self$min_rna & nCount_RNA < self$max_rna)
      return(seurat_object)
    },

    #' Normalize and Stabilize Variance Using SCTransform
    #'
    #' This function normalizes and stabilizes variance in a Seurat object
    #' using the SCTransform method.
    #'
    #' @param seurat_object A Seurat object.
    #' @param verbose       Boolean indicating whether to print detailed
    #'                      messages.
    #'
    #' @return              A normalized Seurat object.
    #' @note                Consider adding `residual_type = "pearson"` to
    #'                      SCTransform if corrected UMI counts are needed.
    normalize_data = function(seurat_object) {
      seurat_object <-
        SCTransform(seurat_object,
                    vars.to.regress = c("nCount_RNA", "percent.mt"),
                    return.only.var.genes = FALSE, verbose = self$verbose,
                    conserve.memory = TRUE)
      return(seurat_object)
    },

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
    annotate_cells_with_singler = function(seurat_object) {
      singler_results <- SingleR(test = seurat_object[["SCT"]]@counts,
                                 ref = self$blueprint_encode,
                                 labels = self$blueprint_encode$label.main)
      seurat_object$blueprint_labels <- singler_results$labels
      seurat_object$blueprint_pvals <- singler_results$scores
      return(seurat_object)
    }
  )
)