#' Compute Distance Matrix
#'
#' Calculates a distance matrix using Pearson correlation.
#'
#' @param dat_mat A matrix of gene expression data (genes x samples).
#'
#' @return        A distance matrix.
compute_distance_matrix <- function(dat_mat) {
  if (!is.matrix(dat_mat)) {
    dat_mat <- as.matrix(dat_mat)
  }

  dist_mat <- as.dist(1 - cor(dat_mat, method = "pearson"))
  return(dist_mat)
}

#' Create Multiple Directories
#'
#' This function creates multiple directories if they do not already exist.
#'
#' @param dir_paths A character vector of directory paths to be created.
#'
#' @return None.
create_directories <- function(dir_paths) {
  for (dir_path in dir_paths) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
  }
}