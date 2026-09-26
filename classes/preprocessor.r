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
    do_filter = TRUE,
    conserve_memory = FALSE,
    variable_features_n = 16000,
    sct_min_cells = 0,
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
    #' @param do_filter    Whether to apply the cell filter. Set FALSE when the
    #'                     input has already been QC'd upstream.
    #'
    #'                     ONLY the filter is skipped. percent.mt is still
    #'                     computed - SCTransform regresses it out, so it is
    #'                     required either way - and SCTransform and SingleR
    #'                     still run, since the integration needs the SCT
    #'                     assay. Deleting this class is therefore not an
    #'                     option; the filter is the only redundant part.
    #'
    #'                     Filtering twice is harmless when the thresholds
    #'                     match and destructive when they do not. Measured on
    #'                     the neutrophil object, which was already filtered at
    #'                     mt<20 / 300<counts<15000 / genes>200: re-filtering at
    #'                     min_rna = 300 drops 0 cells, while this class's old
    #'                     default of 1000 drops 13,101 (25.9%), including >80%
    #'                     of six studies and >50% of sixteen of twenty-five -
    #'                     because sequencing depth here is a dataset property.
    #'                     Defaults to TRUE, so behaviour is unchanged for
    #'                     un-QC'd input.
    #' @param conserve_memory Passed to SCTransform. MUST BE FALSE for this
    #'                     pipeline. Seurat sets return.only.var.genes = TRUE
    #'                     whenever conserve.memory = TRUE, so each object's SCT
    #'                     scale.data keeps only its own ~3,000 variable genes -
    #'                     the explicit return.only.var.genes = FALSE above is
    #'                     silently overridden. PrepSCTIntegration can then only
    #'                     use features present in EVERY object's SCT model, and
    #'                     the intersection of 25 per-object variable-gene sets
    #'                     collapsed to 314 genes on the neutrophil arm (job
    #'                     32455). VIPER got a 314-gene signature holding 2.5% of
    #'                     regulon targets and returned nothing usable.
    #'
    #'                     With it FALSE, residuals are computed for all genes,
    #'                     so the integrated assay spans the full space. Costs
    #'                     roughly n_genes x n_cells x 8 bytes across the list -
    #'                     about 6.3 GB for 16,126 x 49,142. Set TRUE only if
    #'                     that does not fit, and then expect to lose most of
    #'                     the regulon.
    #' @param variable_features_n Passed to SCTransform. Seurat's default of
    #'                     3000 flags only the top 3,000 genes as variable, and
    #'                     SelectIntegrationFeatures draws its anchor candidates
    #'                     from exactly that per-object set - so the pool the
    #'                     anchors come from, and the shared space they must
    #'                     agree on, is capped there. 16000 makes effectively
    #'                     the whole gene space a candidate, which is what this
    #'                     project's own seurat_integration_*.R scripts already
    #'                     do (they use 15000).
    #'
    #'                     This does NOT change what lands in scale.data -
    #'                     return.only.var.genes = FALSE with
    #'                     conserve_memory = FALSE already keeps residuals for
    #'                     every gene. It changes which genes are eligible to
    #'                     become anchors. Defaults to 16000.
    #' @param sct_min_cells Passed to SCTransform as min_cells. Seurat's default
    #'                     of 5 drops genes seen in fewer than 5 cells WITHIN
    #'                     each object, so the shared gene space available to
    #'                     PrepSCTIntegration is the INTERSECTION of 25 (or 85)
    #'                     separately-filtered sets. Measured on the neutrophil
    #'                     Origins, all of which carry the same 16,126-gene axis:
    #'                       min_cells = 5 ->    314 genes   (jobs 32455, 37379)
    #'                                   3 ->    525
    #'                                   1 ->  1,377
    #'                                   0 -> 16,126
    #'                     314 genes left VIPER with 2.5% of its regulon targets
    #'                     and no usable output. Gene filtering belongs upstream
    #'                     where it applies to the whole cohort at once, not
    #'                     per object. Defaults to 0.
    #'
    #'                     Consequence to expect: 14,749 of 16,126 genes have
    #'                     zero counts in at least one neutrophil Origin. Those
    #'                     have no variance to fit, so SCTransform will warn and
    #'                     may return NaN residuals for them in that object.
    #'                     They will not be chosen as anchors - the effective
    #'                     anchor pool sits between the min_cells = 1 figure and
    #'                     the full axis.
    initialize = function(seurat_list, mt_threshold = 25,
                          min_rna = 1000, max_rna = 15000,
                          verbose = FALSE, do_filter = TRUE,
                          conserve_memory = FALSE,
                          variable_features_n = 16000,
                          sct_min_cells = 0) {
      self$seurat_list <- seurat_list
      self$mt_threshold <- mt_threshold
      self$min_rna <- min_rna
      self$max_rna <- max_rna
      self$verbose <- verbose
      self$do_filter <- do_filter
      self$conserve_memory <- conserve_memory
      self$variable_features_n <- variable_features_n
      self$sct_min_cells <- sct_min_cells
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
      if (self$do_filter) {
        n_before <- ncol(seurat_object)
        seurat_object <- private$filter_cells(seurat_object)
        message(sprintf("  filter: %d -> %d cells (%d removed)",
                        n_before, ncol(seurat_object),
                        n_before - ncol(seurat_object)))
      } else {
        message("  filter: skipped (do_filter = FALSE); input already QC'd")
      }
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
                    conserve.memory = self$conserve_memory,
                    variable.features.n = self$variable_features_n,
                    min_cells = self$sct_min_cells)
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