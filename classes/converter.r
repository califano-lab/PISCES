library(R6)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tools)

#' Converter Class for converting Non-Read10X data to Read10X data
#'
#' This R6 class provides methods to convert Non-Read10X data to Read10X data.
#'
#' @field base_data_path The base path to the data directory containing
#'                       patient files and directories.
Converter <- R6Class( # nolint
  "Converter",
  public = list(
    base_data_path = NULL,

    #' Initialize the Converter
    #'
    #' @param base_data_path The base path to the data directory containing
    #'                       patient files and directories.
    initialize = function(base_data_path) {
      self$base_data_path <- base_data_path
    },

    #' Convert ENSEMBL IDs to Gene Names
    #'
    #' This method processes all .rds files in the specified directory,
    #' converts the ENSEMBL IDs to gene names, and saves the processed data.
    #'
    #' @return None.
    convert_ensembl_to_gene_names = function() {
      files <-
        list.files(self$base_data_path, pattern = "\\.rds$", full.names = TRUE)

      for (file in files) {
        tryCatch({
          rds_data <- readRDS(file)
          exprs_data <- rds_data@assayData$exprs
          gene_name_data <- private$ensembl_to_gene_name(exprs_data)

          output_file <- paste0(tools::file_path_sans_ext(file), ".rds")
          saveRDS(gene_name_data, output_file)
        }, error = function(e) {
          message("Failed to process file: ", file, ". Error: ", e$message)
          message(paste0("Please ensure the file is a valid .rds file and",
                         " that it contains the expected data structure."))
          stop(paste0("Also, make sure you are calling the correct conversion",
                      " method in the Converter class."))
        })
      }
    }
  ),

  #############################################################################
  #                           PRIVATE METHODS                                 #
  #############################################################################
  private = list(
    #' Convert ENSEMBL IDs to Gene Names for a Data Frame
    #'
    #' This function converts ENSEMBL IDs to gene names for a given data frame.
    #'
    #' @param dat A data frame with ENSEMBL IDs as row names.
    #'
    #' @return A data frame with gene names as row names.
    ensembl_to_gene_name = function(dat) {
      x <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(dat),
                                 columns = c("SYMBOL"), keytype = "ENSEMBL")
      x <- x[which(!is.na(x[, 2])), ]
      entrezids <- c()
      exclude <- c()
      include <- c()
      for (i in rownames(dat)) {
        ids <- x[which(x[, 1] == i), 2]
        if (length(ids) == 0) {
          exclude <- c(exclude, i)
        } else {
          entrezids <- c(entrezids, ids[1])
          include <- c(include, i)
        }
      }
      dat <- dat[include, ]
      rownames(dat) <- make.unique(entrezids)
      return(dat)
    }
  )
)
