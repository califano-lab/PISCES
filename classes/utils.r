Utils <- R6Class( # nolint
  "Utils",
  public = list(

    init = function(base_output_path, plot_output_path, aracne_output_path,
                    viper_output_path, base_data_path, aracne_binary_path,
                    regulator_dir_path) {
      plot_output_path <- file.path(base_output_path, "plots")
      aracne_output_path <- file.path(base_output_path, "aracne_results")
      viper_output_path <- file.path(base_output_path, "viper_results")

      self$create_directories(c(
        base_output_path,
        plot_output_path,
        aracne_output_path,
        viper_output_path
      ))

      self$do_paths_exist(
        c(base_data_path, base_output_path,
          aracne_binary_path, regulator_dir_path)
      )

      return(list(
        plot_output_path = plot_output_path,
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
  )
)