source("classes/utils.r")

library(Seurat)

#' Runner Class for Running ARACNe and VIPER Analyses
#'
#' This class provides methods to run ARACNe and VIPER analyses on a given
#' metacell matrix and Seurat object. It initializes the necessary parameters
#' and includes methods for generating regulon objects from ARACNe output
#' and running VIPER analysis.
#'
#' @field metacell_mat The metacell matrix.
#' @field base_output_path The base path for output files.
#' @field aracne_binary_path The path to the ARACNe binary executable.
#' @field aracne_output_path The path to save ARACNe output files.
#' @field regulator_files A named list of paths to regulator files.
#' @field seurat_obj A Seurat object containing the expression matrix.
#' @field viper_output_path The path to save VIPER output files.
#' @field threads Number of threads to use for ARACNe computation.
#' @field seed Seed for random number generation to ensure reproducibility.
Runner <- R6Class( # nolint
  "Runner",
  public = list(
    metacell_mat = NULL,
    base_output_path = NULL,
    aracne_binary_path = NULL,
    aracne_output_path = NULL,
    regulator_files = NULL,
    seurat_obj = NULL,
    viper_output_path = NULL,
    threads = NULL,
    seed = NULL,
    utils = Utils$new(),

    #' Initialize the Runner Class
    #'
    #' @param metacell_mat       The metacell matrix.
    #' @param base_output_path   The base path for output files.
    #' @param aracne_binary_path The path to the ARACNe binary executable.
    #' @param aracne_output_path The path to save ARACNe output files.
    #' @param regulator_files    A named list of paths to regulator files.
    #' @param seurat_obj         A Seurat object containing the expression
    #'                           matrix.
    #' @param viper_output_path  The path to save VIPER output files.
    #' @param threads            Number of threads to use for ARACNe
    #'                           computation.
    #' @param seed               Seed for random number generation to ensure
    #'                           reproducibility.
    initialize = function(metacell_mat, base_output_path, aracne_binary_path,
                          aracne_output_path, regulator_files, seurat_obj,
                          viper_output_path, threads = 4, seed = 42) {
      self$metacell_mat <- metacell_mat
      self$base_output_path <- base_output_path
      self$aracne_binary_path <- aracne_binary_path
      self$aracne_output_path <- aracne_output_path
      self$regulator_files <- regulator_files
      self$seurat_obj <- seurat_obj
      self$viper_output_path <- viper_output_path
      self$threads <- threads
      self$seed <- seed
    },

    #' Run ARACNe Analysis
    #'
    #' This function runs the ARACNe analysis on the metacell matrix.
    run_aracne = function() {
      if (length(self$metacell_mat) <= 0) {
        stop("Metacell matrix is empty.")
      }
      expression_files <-
        private$prep_and_save_expr_for_aracne(self$base_output_path,
                                              "_all_all.txt.tsv")
      private$run_aracne_helper(self$aracne_binary_path, expression_files,
                                self$regulator_files, self$aracne_output_path,
                                threads = self$threads, seed = self$seed)
    },

    #' Run VIPER Analysis
    #'
    #' This function runs the VIPER analysis using the generated regulon
    #' objects from the ARACNe output.
    run_viper = function() {
      exp_mat <- private$load_expr_matrix(self$seurat_obj)
      regulon_list <-
        private$generate_regulon_objects(self$aracne_output_path, exp_mat,
                                         self$base_output_path)
      viper_results <- private$run_viper_helper(exp_mat, regulon_list)
      viper_results_path <-
        file.path(self$viper_output_path, "viper_results.rds")
      save_viper_results(viper_results, viper_results_path)
      message("VIPER analysis completed and results saved.")
    }
  ),

  private = list(
    #' Prepare and Save Expression Matrix for ARACNe Analysis
    #'
    #' This function prepares and saves an expression matrix for use in ARACNe
    #' analysis. It first removes duplicate genes, then subsets the matrix if
    #' necessary, and finally saves the matrix in both RDS and TSV formats.
    #'
    #' @param base_output_path  The base path for the output files.
    #' @param file_suffix       The suffix to append to the input file name.
    prep_and_save_expr_for_aracne = function(base_output_path, file_suffix) {
      files <- list.files(base_output_path, pattern = paste0(file_suffix, "$"),
                          full.names = TRUE)
      expr_files <- lapply(files, function(file_path) {
        expr_data <- read.table(file_path, header = TRUE, sep = "\t",
                                row.names = 1)
        out_path <- gsub(".tsv$", "_for_aracne.tsv", file_path)
        private$save_matrix_for_aracne(expr_data, out_path)
        return(out_path)
      })
      return(expr_files)
    },

    #' Save Expression Matrix for ARACNe
    #'
    #' Saves an expression matrix in a format required by ARACNe3.
    #'
    #' @param expression_matrix The normalized expression matrix to be saved.
    #' @param output_file       Path for the output file.
    save_matrix_for_aracne = function(expression_matrix, output_file) {
      # Ensure the matrix has proper row and column names
      if (is.null(colnames(expression_matrix)) ||
            is.null(rownames(expression_matrix))) {
        stop("Expression matrix must have row and column names.")
      }

      # Ensure there is an extra column in the header (if not already present)
      cat("\t", file = output_file, append = FALSE)
      suppressWarnings(
        write.table(expression_matrix, file = output_file, sep = "\t",
                    quote = FALSE, row.names = TRUE, col.names = TRUE,
                    append = TRUE)
      )
    },

    #' Run ARACNe for Each Expression and Regulator File
    #'
    #' This function runs ARACNe for each combination of expression and
    #' regulator files, creating the necessary output directories and executing
    #' the ARACNe command.
    #'
    #' @param aracne_bin       Path to the ARACNe binary executable.
    #' @param expression_files A list of paths to expression files.
    #' @param regulator_files  A named list of paths to regulator files.
    #' @param output_base_dir  The base directory where ARACNe output will be
    #'                         saved.
    #' @param threads          Number of threads to use for ARACNe.
    #' @param seed             Seed for random number generation to ensure
    #'                         reproducibility.
    #'
    #' @return                 None. The function executes ARACNe and saves the
    #'                         results to the specified output directories.
    run_aracne_helper = function(aracne_bin, expression_files, regulator_files,
                                 output_base_dir, threads, seed) {
      for (reg_name in names(regulator_files)) {
        regulator_file <- regulator_files[[reg_name]]
        for (exp_file in expression_files) {
          exp_file_base <-
            gsub("_all_all.txt_for_aracne.tsv", "", basename(exp_file))
          output_dir <-
            file.path(output_base_dir, paste0(reg_name, "_", exp_file_base))
          self$utils$create_directories(list(output_dir))
          private$execute_aracne(aracne_bin, exp_file, regulator_file,
                                 output_dir, threads, seed)
        }
      }
    },

    #' Executes ARACNe3 on the provided expression matrix and regulator list.
    #'
    #' @param aracne_bin      Path to the ARACNe3 binary.
    #' @param exp_file        Path to the expression matrix file.
    #' @param regulators_file Path to the file containing regulator gene names.
    #' @param output_dir      Directory to store ARACNe output.
    #' @param threads         Number of threads to use for ARACNe computation.
    #' @param seed            Seed for random number generation in ARACNe.
    execute_aracne = function(aracne_bin, exp_file, regulators_file,
                              output_dir, threads, seed) {
      cmd <-
        sprintf("%s -e %s -r %s -o %s --threads %d --seed %d",
                aracne_bin, exp_file, regulators_file, output_dir, threads,
                seed)
      system(cmd)
    },

    #' Load Expression Matrix from Seurat Object
    #'
    #' This function loads the expression matrix from a Seurat object.
    #'
    #' @param seurat_obj A Seurat object containing the expression matrix.
    #'
    #' @return           The expression matrix from the Seurat object.
    load_expr_matrix = function(seurat_obj) {
      exp_mat <-
        GetAssayData(object = seurat_obj, assay = "integrated",
                     layer = "scale.data")
      if (is.null(exp_mat) || ncol(exp_mat) == 0 || nrow(exp_mat) == 0) {
        stop(paste0("Expression matrix is empty or NULL.",
                    " Check your Seurat object and data extraction steps."))
      }
      exp_mat <- as.matrix(exp_mat)
      return(exp_mat)
    },

    #' Process ARACNe Output Files to Generate Regulon Objects
    #'
    #' This function processes ARACNe output files to generate regulon objects
    #' suitable for VIPER analysis. It loads ARACNe output files, prepares the
    #' data, and generates regulon objects.
    #'
    #' @param aracne_output_base_dir The base directory containing ARACNe
    #'                               output files.
    #' @param exp_mat                The expression matrix used to generate the
    #'                               ARACNe network (genes x samples).
    #' @param output_base_path       The base directory where the processed
    #'                               regulon objects will be saved.
    #'
    #' @return                       A list of regulon objects generated from
    #'                               the ARACNe output files.
    generate_regulon_objects = function(aracne_output_base_dir, exp_mat,
                                        output_base_path) {
      aracne_output_files <-
        private$get_all_aracne_files(aracne_output_base_dir)

      if (length(aracne_output_files) == 0) {
        stop("No ARACNe output files found in the directory.")
      }

      regulon_list <- lapply(aracne_output_files, function(aracne_file) {
        aracne_data_for_viper <-
          private$prep_aracne_output_for_viper(aracne_file)
        prefix <- gsub("-metaCells$", "", basename(dirname(aracne_file)))
        regulon <- private$generate_regulon(aracne_data_for_viper, exp_mat,
                                            output_base_path, prefix)
        return(regulon)
      })

      return(regulon_list)
    },

    #' Aggregate ARACNe Output Files from All Directories
    #'
    #' This helper function aggregates ARACNe output files from all directories
    #' within the specified base directory, matching the provided file name
    #' pattern.
    #'
    #' @param base_dir The base directory containing the ARACNe output files.
    #' @param pattern  A regular expression pattern to match ARACNe output
    #'                 files.
    #'
    #' @return         A character vector of file paths to the ARACNe output
    #'                 files.
    get_all_aracne_files =
      function(base_dir, pattern = "consolidated-net_.*\\.tsv$") {
        all_files <- list.files(base_dir, pattern = pattern, full.names = TRUE,
                                recursive = TRUE)
        return(all_files)
      },

    #' Load and Process ARACNe Output File
    #'
    #' This helper function loads an ARACNe output file, processes its content,
    #' and prepares it for VIPER analysis.
    #'
    #' @param aracne_file The path to the ARACNe output file.
    #'
    #' @return            A data frame containing the processed ARACNe data with
    #'                    columns "regulator", "target", and "mi"
    #'                    (mutual information).
    prep_aracne_output_for_viper = function(aracne_file) {
      cat("Loading ARACNe output file:", aracne_file, "\n")

      # Load ARACNe output file without headers
      aracne_data <- read.table(aracne_file, header = FALSE, sep = "\t",
                                check.names = FALSE, stringsAsFactors = FALSE,
                                skip = 1)

      # Only include the first three columns
      aracne_data <- aracne_data[, 1:3]

      # Define column names manually
      colnames(aracne_data) <- c("regulator", "target", "mi")

      # Convert the 'mi' column to numeric
      aracne_data$mi <- as.numeric(aracne_data$mi)

      return(aracne_data)
    },

    #' Generate Regulon Object from ARACNe Output
    #'
    #' This helper function generates a regulon object from ARACNe output data
    #' and an expression matrix. The regulon object is then saved to the
    #' specified output directory.
    #'
    #' @param aracne_data      A data frame containing ARACNe output data.
    #' @param exp_mat          An expression matrix associated with the ARACNe
    #'                         network.
    #' @param output_base_path The base directory where the regulon objects
    #'                         will be saved.
    #' @param file_prefix      The prefix for the saved regulon files.
    #'
    #' @return                 A pruned regulon object suitable for VIPER
    #'                         analysis.
    generate_regulon = function(aracne_data, exp_mat, output_base_path,
                                file_prefix) {
      # Process ARACNe results for VIPER analysis
      private$reg_process(aracne_data, exp_mat, output_base_path, file_prefix)

      # Load the pruned regulon object
      pruned_regulon_file <-
        file.path(output_base_path, paste0(file_prefix, "_pruned.rds"))
      pruned_regulon <- readRDS(pruned_regulon_file)

      return(pruned_regulon)
    },

    #' Run VIPER on a List of Regulon Objects
    #'
    #' This function runs VIPER analysis on a list of regulon objects using an
    #' expression matrix. It processes each regulon object and returns the
    #' VIPER scores.
    #'
    #' @param exp_mat      An expression matrix used for VIPER analysis.
    #' @param regulon_list A list of regulon objects generated from ARACNe
    #'                     output.
    #'
    #' @return             A list of VIPER results for each regulon object.
    run_viper_helper = function(exp_mat, regulon_list) {
      viper_results <- lapply(regulon_list, function(regulon) {
        viper_scores <- private$execute_viper(exp_mat, regulon)
        return(viper_scores)
      })

      # Check contents of viper_results
      if (length(viper_results) == 0 || any(sapply(viper_results, is.null))) {
        stop("VIPER results are empty or not properly formed.")
      }

      return(viper_results)
    },

    #' Process ARACNe Results for VIPER Analysis
    #'
    #' Converts ARACNe output into a regulon object suitable for VIPER analysis,
    #' including an optional pruning step to refine the regulon.
    #'
    #' @param a_file   Path to the ARACNe final network file in TSV format.
    #' @param exp_mat  Expression matrix used to generate the ARACNe network
    #'                 (genes x samples).
    #' @param out_dir  Directory where the processed regulon objects will be
    #'                 saved.
    #' @param out_name Prefix for the saved regulon files.
    reg_process = function(a_file, exp_mat, out_dir, out_name = "") {
      require(viper)

      # Convert ARACNe output to regulon object
      processed_reg <- private$convert_to_regulon(a_file, exp_mat)

      # Save the unpruned regulon object
      private$save_regulon(processed_reg, out_dir,
                           paste0(out_name, "_unpruned.rds"))

      # Prune the regulon to refine it
      pruned_reg <- private$prune_regulon(processed_reg)

      # Save the pruned regulon object
      private$save_regulon(pruned_reg, out_dir,
                           paste0(out_name, "_pruned.rds"))
    },

    #' Convert ARACNe Output to Regulon Object
    #'
    #' @param a_file  Path to the ARACNe network file.
    #' @param exp_mat Expression matrix associated with the ARACNe network.
    #'
    #' @return        A regulon object suitable for VIPER analysis.
    convert_to_regulon = function(aracne_data, exp_mat) {
      if (is.null(aracne_data) || nrow(aracne_data) == 0) {
        stop("ARACNe data is empty or not available.")
      }

      if (is.null(dim(exp_mat))) {
        stop("Expression matrix is not correctly formatted or is NULL.")
      }

      # Create a temporary file to store processed ARACNe data
      temp_file <- tempfile()
      write.table(aracne_data, temp_file, sep = "\t", row.names = FALSE,
                  col.names = FALSE, quote = FALSE)

      tryCatch({
        regulon_object <- aracne2regulon(afile = temp_file, eset = exp_mat,
                                         format = "3col", verbose = TRUE)
      }, error = function(e) {
        cat("Error during regulon conversion: ", e$message, "\n")
        stop("Failed to convert ARACNe output to regulon object: ", e$message)
      })

      unlink(temp_file)

      if (is.null(regulon_object) || length(regulon_object) == 0) {
        stop("Regulon object is NULL or empty.")
      }

      return(regulon_object)
    },

    #' Save Regulon Object to File
    #'
    #' @param regulon   The regulon object to be saved.
    #' @param out_dir   The directory for saving the regulon file.
    #' @param file_name The name of the file to save the regulon object in.
    save_regulon = function(regulon, out_dir, file_name) {
      saveRDS(regulon, file = file.path(out_dir, file_name))
    },

    #' Prune Regulon Object
    #'
    #' Applies pruning to a regulon object to refine its content.
    #'
    #' @param regulon The regulon object to be pruned.
    #'
    #' @return        A pruned regulon object.
    prune_regulon = function(regulon) {
      viper::pruneRegulon(regulon, 50, adaptive = FALSE, eliminate = TRUE)
    },

    #' Run VIPER Analysis on a Single Regulon
    #'
    #' This function executes VIPER analysis on a single regulon using an
    #' expression matrix. It handles errors during the execution and returns the
    #' VIPER scores.
    #'
    #' @param exp_mat An expression matrix used for VIPER analysis.
    #' @param regulon A regulon object generated from ARACNe output.
    #'
    #' @return VIPER scores for the given regulon object or NULL if an error
    #'         occurs.
    execute_viper = function(exp_mat, regulon) {
      viper_scores <- tryCatch({
        viper(exp_mat, regulon)
      }, error = function(e) {
        cat("Error in VIPER analysis:", e$message, "\n")
        NULL
      })

      return(viper_scores)
    },

    #' Save VIPER Results to File
    #'
    #' This function saves the VIPER analysis results to a specified file in
    #' RDS format.
    #'
    #' @param viper_results A list of VIPER results to be saved.
    #' @param output_path   The file path where the VIPER results will be
    #'                      saved.
    #'
    #' @return              None.
    save_viper_results = function(viper_results, output_path) {
      saveRDS(viper_results, file = output_path)
    }
  )
)