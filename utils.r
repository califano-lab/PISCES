Utils <- R6Class( # nolint
  "Utils",
  public = list(
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

    #' Notify Non-Existent Directories
    #'
    #' This function checks a list of directory paths and notifies the user
    #' which paths do not exist.
    #'
    #' @param dir_paths A character vector of directory paths to be checked.
    #'
    #' @return A character vector of non-existent directory paths.
    do_directories_exist = function(dir_paths) {
      non_existent_dirs <- dir_paths[!dir.exists(dir_paths)]
      if (length(non_existent_dirs) > 0) {
        message("The following directory paths do not exist:")
        print(non_existent_dirs)
      } else {
        message("All directory paths exist.")
      }
      return(non_existent_dirs)
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