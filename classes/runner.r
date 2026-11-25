source("classes/utils.r")

library(Seurat)

#' Runner Class for Running ARACNe and VIPER Analyses
#'
#' This class runs ARACNe3 on the metacell matrices produced by
#' `MetacellGenerator$generate_metacell_matrices()`, converts the resulting
#' networks into regulons, and runs VIPER with them.
#'
#' ARACNe3 is run once per (regulator set, metacell matrix) pair, on exactly
#' the matrices passed to `initialize()`: the subsampled, CPM-normalised
#' metacells of each cluster. `run_viper()` reads only the networks of those
#' same pairs, so outputs left in the directories by earlier runs are never
#' mixed in.
#'
#' @field metacell_mat Named list of metacell matrices (genes x metacells),
#'        as returned by `MetacellGenerator$generate_metacell_matrices()`.
#' @field base_output_path The base path for output files.
#' @field aracne_binary_path The path to the ARACNe binary executable.
#' @field aracne_output_path The path to save ARACNe output files.
#' @field regulator_files A named list of paths to regulator files.
#' @field seurat_obj A Seurat object containing the expression matrix.
#' @field viper_output_path The path to save VIPER output files.
#' @field threads Number of threads to use for ARACNe computation.
#' @field seed Seed for random number generation to ensure reproducibility.
#' @field num_subnets Number of ARACNe3 subnetworks per network (ARACNe3's
#'        `-x` option).
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
    num_subnets = NULL,
    utils = Utils$new(),

    #' Initialize the Runner Class
    #'
    #' @param metacell_mat       Named list of metacell matrices, as returned
    #'                           by the MetacellGenerator method
    #'                           `generate_metacell_matrices()`. Each matrix is
    #'                           used as-is as one ARACNe input.
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
    #' @param num_subnets        Number of subnetworks ARACNe3 builds and
    #'                           consolidates into each network (its `-x`
    #'                           option). Each subnetwork is inferred from a
    #'                           different random 63% of the metacells. With 1
    #'                           (ARACNe3's default) every edge of the
    #'                           consolidated network has the same p-value;
    #'                           with more, edges are scored by how many
    #'                           subnetworks recover them, at a proportional
    #'                           cost in run time. VIPER uses only the mutual
    #'                           information column, so it runs either way.
    initialize = function(metacell_mat, base_output_path, aracne_binary_path,
                          aracne_output_path, regulator_files, seurat_obj,
                          viper_output_path, threads = 4, seed = 42,
                          num_subnets = 1) {
      self$utils$check_whole_number(threads, "threads")
      self$utils$check_whole_number(seed, "seed", minimum = 0)
      self$utils$check_whole_number(num_subnets, "num_subnets")

      self$metacell_mat <- metacell_mat
      self$base_output_path <- base_output_path
      self$aracne_binary_path <- aracne_binary_path
      self$aracne_output_path <- aracne_output_path
      self$regulator_files <- regulator_files
      self$seurat_obj <- seurat_obj
      self$viper_output_path <- viper_output_path
      self$threads <- threads
      self$seed <- seed
      self$num_subnets <- num_subnets
    },

    #' Run ARACNe Analysis
    #'
    #' Writes each metacell matrix to
    #' `{aracne_output_path}/{name}_for_aracne.tsv` and runs ARACNe3 on it once
    #' per regulator set. The network for regulator set `reg` and matrix
    #' `name` is written to `{aracne_output_path}/{reg}_{name}/`, which is
    #' cleared first so that no files from an earlier run remain in it.
    run_aracne = function() {
      private$check_metacell_mats(self$metacell_mat)
      private$check_regulator_files(self$regulator_files)
      self$utils$create_directories(self$aracne_output_path)

      for (mat_name in names(self$metacell_mat)) {
        exp_file <- private$aracne_input_file(mat_name)
        private$save_matrix_for_aracne(self$metacell_mat[[mat_name]],
                                       exp_file)

        for (reg_name in names(self$regulator_files)) {
          output_dir <- private$aracne_run_dir(reg_name, mat_name)
          unlink(output_dir, recursive = TRUE)
          self$utils$create_directories(output_dir)
          private$execute_aracne(self$aracne_binary_path, exp_file,
                                 self$regulator_files[[reg_name]], output_dir,
                                 self$threads, self$seed, self$num_subnets)
        }
      }
    },

    #' Run VIPER Analysis
    #'
    #' Converts the ARACNe networks produced by `run_aracne()` for this
    #' Runner's metacell matrices and regulator sets into pruned regulons
    #' (saved to `base_output_path`), runs VIPER with all of them on the
    #' Seurat object's integrated scale.data, and saves the result as
    #' `{viper_output_path}/viper_results.rds`.
    run_viper = function() {
      private$check_metacell_mats(self$metacell_mat)
      private$check_regulator_files(self$regulator_files)

      exp_mat <- private$load_expr_matrix(self$seurat_obj)
      regulon_list <-
        private$generate_regulon_objects(private$expected_network_files(),
                                         exp_mat, self$base_output_path)
      viper_results <- private$run_viper_helper(exp_mat, regulon_list)
      viper_results_path <-
        file.path(self$viper_output_path, "viper_results.rds")
      private$save_viper_results(viper_results, viper_results_path)
      message("VIPER analysis completed and results saved.")
    }
  ),

  #############################################################################
  #                           PRIVATE METHODS                                 #
  #############################################################################
  private = list(
    #' Check the Metacell Matrices Before Using Them
    #'
    #' Each element must be a numeric genes x metacells matrix with unique
    #' gene names (ARACNe3 rejects duplicates), and the list must be named
    #' uniquely, because the names identify each matrix's files.
    #'
    #' @param metacell_mats The value of the `metacell_mat` field.
    check_metacell_mats = function(metacell_mats) {
      usage <- paste("`metacell_mat` must be the named list returned by",
                     "MetacellGenerator$generate_metacell_matrices().")
      mat_names <- names(metacell_mats)
      if (!is.list(metacell_mats) || length(metacell_mats) == 0) {
        stop("Metacell matrix list is empty. ", usage)
      }
      if (is.null(mat_names) || any(!nzchar(mat_names)) ||
            anyDuplicated(mat_names)) {
        stop("Every metacell matrix needs a unique name. ", usage)
      }
      for (mat_name in mat_names) {
        mat <- metacell_mats[[mat_name]]
        if (!is.matrix(mat) || !is.numeric(mat) || is.null(rownames(mat)) ||
              is.null(colnames(mat))) {
          stop(sprintf(paste("Metacell matrix '%s' must be a numeric matrix",
                             "with gene (row) and metacell (column) names."),
                       mat_name))
        }
        if (anyDuplicated(rownames(mat))) {
          stop(sprintf("Metacell matrix '%s' has duplicate gene names.",
                       mat_name))
        }
      }
    },

    #' Check the Regulator Files
    #'
    #' @param regulator_files The value of the `regulator_files` field.
    check_regulator_files = function(regulator_files) {
      reg_names <- names(regulator_files)
      if (length(regulator_files) == 0 || is.null(reg_names) ||
            any(!nzchar(reg_names)) || anyDuplicated(reg_names)) {
        stop("`regulator_files` must be a non-empty list with unique names.")
      }
      missing <- unlist(regulator_files)[!file.exists(unlist(regulator_files))]
      if (length(missing) > 0) {
        stop("Regulator file(s) not found: ", paste(missing, collapse = ", "))
      }
    },

    #' Path of the ARACNe3 Input File Written for a Metacell Matrix
    #'
    #' @param mat_name Name of the metacell matrix.
    aracne_input_file = function(mat_name) {
      file.path(self$aracne_output_path, paste0(mat_name, "_for_aracne.tsv"))
    },

    #' ARACNe3 Output Directory for a Regulator Set and Metacell Matrix
    #'
    #' The directory name also prefixes the regulon files written by
    #' `run_viper()`.
    #'
    #' @param reg_name Name of the regulator set.
    #' @param mat_name Name of the metacell matrix.
    aracne_run_dir = function(reg_name, mat_name) {
      file.path(self$aracne_output_path, paste0(reg_name, "_", mat_name))
    },

    #' Consolidated Network File Inside an ARACNe3 Output Directory
    #'
    #' ARACNe3 names the file after its run id, which Runner leaves at the
    #' ARACNe3 default, "defaultid".
    #'
    #' @param output_dir An ARACNe3 output directory.
    consolidated_network_file = function(output_dir) {
      file.path(output_dir, "consolidated-net_defaultid.tsv")
    },

    #' Network Files This Runner's ARACNe Runs Produce
    #'
    #' @return A character vector of consolidated network paths, named by the
    #'         output directory name (`{reg}_{name}`), for every regulator set
    #'         and metacell matrix.
    expected_network_files = function() {
      run_dirs <- unlist(lapply(names(self$regulator_files), function(reg) {
        vapply(names(self$metacell_mat), function(mat_name) {
          private$aracne_run_dir(reg, mat_name)
        }, character(1), USE.NAMES = FALSE)
      }))
      network_files <- vapply(run_dirs, private$consolidated_network_file,
                              character(1))
      names(network_files) <- basename(run_dirs)
      return(network_files)
    },

    #' Save Expression Matrix for ARACNe
    #'
    #' Writes a tab-separated (genes + 1) x (samples + 1) file: a header line
    #' with an empty first cell followed by the sample names, then one line
    #' per gene. ARACNe3 counts the samples from the tabs in the header line,
    #' so the empty first cell (the leading tab) is required.
    #'
    #' @param expression_matrix The normalized expression matrix to be saved.
    #' @param output_file       Path for the output file.
    save_matrix_for_aracne = function(expression_matrix, output_file) {
      # Ensure the matrix has proper row and column names
      if (is.null(colnames(expression_matrix)) ||
            is.null(rownames(expression_matrix))) {
        stop("Expression matrix must have row and column names.")
      }

      # write.table() omits the corner cell when writing row names, so write
      # it (an empty field) first.
      cat("\t", file = output_file, append = FALSE)
      suppressWarnings(
        write.table(expression_matrix, file = output_file, sep = "\t",
                    quote = FALSE, row.names = TRUE, col.names = TRUE,
                    append = TRUE)
      )
    },

    #' Executes ARACNe3 on the provided expression matrix and regulator list.
    #'
    #' Stops with an error if ARACNe3 fails or does not write a consolidated
    #' network, rather than continuing to VIPER with a missing network.
    #'
    #' @param aracne_bin      Path to the ARACNe3 binary.
    #' @param exp_file        Path to the expression matrix file.
    #' @param regulators_file Path to the file containing regulator gene names.
    #' @param output_dir      Directory to store ARACNe output.
    #' @param threads         Number of threads to use for ARACNe computation.
    #' @param seed            Seed for random number generation in ARACNe.
    #' @param num_subnets     Number of subnetworks to build (ARACNe3's `-x`).
    execute_aracne = function(aracne_bin, exp_file, regulators_file,
                              output_dir, threads, seed, num_subnets) {
      args <- c("-e", shQuote(exp_file),
                "-r", shQuote(regulators_file),
                "-o", shQuote(output_dir),
                "--threads", format(threads, scientific = FALSE),
                "--seed", format(seed, scientific = FALSE),
                "-x", format(num_subnets, scientific = FALSE))
      status <- suppressWarnings(system2(aracne_bin, args))
      if (!identical(as.integer(status), 0L)) {
        stop(sprintf("ARACNe3 failed (exit status %s): %s %s",
                     status, aracne_bin, paste(args, collapse = " ")))
      }
      if (!file.exists(private$consolidated_network_file(output_dir))) {
        stop(sprintf("ARACNe3 did not write a consolidated network in %s.",
                     output_dir))
      }
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
    #' Converts each of this Runner's ARACNe networks into a pruned regulon,
    #' saved as `{output_base_path}/{reg}_{name}_unpruned.rds` and
    #' `..._pruned.rds`.
    #'
    #' @param network_files    Consolidated network paths named by their
    #'                         output directory (see
    #'                         `expected_network_files()`).
    #' @param exp_mat          The expression matrix used to generate the
    #'                         ARACNe network (genes x samples).
    #' @param output_base_path The base directory where the processed
    #'                         regulon objects will be saved.
    #'
    #' @return                 A list of pruned regulon objects, named like
    #'                         `network_files`.
    generate_regulon_objects = function(network_files, exp_mat,
                                        output_base_path) {
      missing <- network_files[!file.exists(network_files)]
      if (length(missing) > 0) {
        stop("ARACNe network(s) not found; run run_aracne() first:\n",
             paste(missing, collapse = "\n"))
      }

      regulon_list <- lapply(names(network_files), function(prefix) {
        aracne_data_for_viper <-
          private$prep_aracne_output_for_viper(network_files[[prefix]])
        private$generate_regulon(aracne_data_for_viper, exp_mat,
                                 output_base_path, prefix)
      })
      names(regulon_list) <- names(network_files)

      return(regulon_list)
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
      viper_results <- private$execute_viper(exp_mat, regulon_list)

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
      viper::pruneRegulon(regulon, 50, adaptive = TRUE, eliminate = TRUE)
    },

    #' Run VIPER Analysis on a Single Regulon
    #'
    #' This function executes VIPER analysis on a single regulon using an
    #' expression matrix. It handles errors during the execution and returns
    #' the VIPER scores.
    #'
    #' @param exp_mat An expression matrix used for VIPER analysis.
    #' @param regulon_list A regulon list generated from ARACNe output.
    #'
    #' @return VIPER scores for the given regulon object or NULL if an error
    #'         occurs.
    execute_viper = function(exp_mat, regulon_list) {
      viper_scores <- tryCatch({
        viper(exp_mat, regulon_list)
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