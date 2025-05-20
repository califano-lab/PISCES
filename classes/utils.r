library(R6)

#' Utils Class
#'
#' The `Utils` class provides various utility functions for path
#' creation/checking, distance matrix calculation, and other small helper tasks
#' used throughout the pipeline.
Utils <- R6Class( # nolint
  "Utils",
  public = list(

    #' Initialize Utility Paths
    #'
    #' This method initializes paths for plotting, ARACNe, and VIPER output
    #' directories based on the specified base paths. It also checks whether
    #' certain required directories/files exist.
    #'
    #' @param base_data_path     A character string specifying the base path
    #'                           containing patient data or other input files
    #'                           for analysis.
    #' @param base_output_path   A character string specifying the base path
    #'                           for all output files/directories.
    #' @param aracne_binary_path A character string specifying the path to the
    #'                           ARACNe3 binary executable.
    #' @param regulator_dir_path A character string specifying the path to the
    #'                           directory containing regulator .txt files.
    #'
    #' @return A named list containing the following elements:
    #'         - plot_output_path
    #'         - plot_viper_output_path
    #'         - aracne_output_path
    #'         - viper_output_path
    init = function(base_data_path, base_output_path, aracne_binary_path,
                    regulator_dir_path) {
      plot_output_path <- file.path(base_output_path, "plots")
      plot_viper_output_path <-
        file.path(plot_output_path, "reclustered_viper_plots")
      aracne_output_path <- file.path(base_output_path, "aracne_results")
      viper_output_path <- file.path(base_output_path, "viper_results")

      self$create_directories(c(
        base_output_path,
        plot_output_path,
        plot_viper_output_path,
        aracne_output_path,
        viper_output_path
      ))

      self$do_paths_exist(
        c(base_data_path, base_output_path,
          aracne_binary_path, regulator_dir_path)
      )

      return(list(
        plot_output_path = plot_output_path,
        plot_viper_output_path = plot_viper_output_path,
        aracne_output_path = aracne_output_path,
        viper_output_path = viper_output_path
      ))
    },

    #' Compute Distance Matrix
    #'
    #' Calculates a distance matrix using Pearson correlation.
    #'
    #' @param dat_mat A matrix of gene expression data (genes x samples).
    #'
    #' @return        A distance matrix.
    compute_distance_matrix = function(dat_mat) {
      if (!is.matrix(dat_mat)) {
        dat_mat <- as.matrix(dat_mat)
      }

      dist_mat <- as.dist(1 - cor(dat_mat, method = "pearson"))
      return(dist_mat)
    },

    #' Notify Non-Existent Paths
    #'
    #' This function checks a list of paths and notifies the user
    #' which paths do not exist.
    #'
    #' @param paths A character vector of paths to be checked.
    #'
    #' @return      A character vector of non-existent paths.
    do_paths_exist = function(paths) {
      normalized_paths <-
        sapply(paths, normalizePath, winslash = "/", mustWork = FALSE)
      non_existent_paths <- normalized_paths[!file.exists(normalized_paths)]

      if (length(non_existent_paths) > 0) {
        message("The following paths do not exist:")
        stop(non_existent_paths)
      }
      return(non_existent_paths)
    },

    #' Create Multiple Directories
    #'
    #' This function creates multiple directories if they do not already exist.
    #'
    #' @param dir_paths A character vector of directory paths to be created.
    #'
    #' @return None.
    create_directories = function(dir_paths) {
      for (dir_path in dir_paths) {
        if (!dir.exists(dir_path)) {
          dir.create(dir_path, recursive = TRUE)
        }
      }
    }
  ),

  #############################################################################
  #                           PRIVATE METHODS                                 #
  #############################################################################
  private = list()
)