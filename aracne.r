source("utils.r")

#' Prepare and Save Expression Matrix for ARACNe Analysis
#'
#' This function prepares and saves an expression matrix for use in ARACNe
#' analysis. It first removes duplicate genes, then subsets the matrix if
#' necessary, and finally saves the matrix in both RDS and TSV formats.
#'
#' @param base_output_path  The base path for the output files.
#' @param file_suffix       The suffix to append to the input file name.
#'
prep_and_save_expr_for_aracne <- function(base_output_path, file_suffix) {
  files <- list.files(base_output_path, pattern = paste0(file_suffix, "$"),
                      full.names = TRUE)
  expr_files <- lapply(files, function(file_path) {
    expr_data <- read.table(file_path, header = TRUE, sep = "\t",
                            row.names = 1)
    out_path <- gsub(".tsv$", "_for_aracne.tsv", file_path)
    save_matrix_for_aracne(expr_data, out_path)
    return(out_path)
  })
  return(expr_files)
}

#' Save Expression Matrix for ARACNe
#'
#' Saves an expression matrix in a format required by ARACNe3.
#'
#' @param expression_matrix The normalized expression matrix to be saved.
#' @param output_file Path for the output file.
save_matrix_for_aracne <- function(expression_matrix, output_file) {
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
}

#' Run ARACNe for Each Expression and Regulator File
#'
#' This function runs ARACNe for each combination of expression and regulator
#' files, creating the necessary output directories and executing the ARACNe
#' command.
#'
#' @param aracne_bin Path to the ARACNe binary executable.
#' @param expression_files A list of paths to expression files.
#' @param regulator_files A named list of paths to regulator files.
#' @param output_base_dir The base directory where ARACNe output will be saved.
#' @param threads Number of threads to use for ARACNe.
#' @param seed Seed for random number generation to ensure reproducibility.
#'
#' @return None. The function executes ARACNe and saves the results to the
#'         specified output directories.
run_aracne <- function(aracne_bin, expression_files, regulator_files,
                       output_base_dir, threads, seed) {
  for (reg_name in names(regulator_files)) {
    regulator_file <- regulator_files[[reg_name]]
    for (exp_file in expression_files) {
      exp_file_base <-
        gsub("_all_all.txt_for_aracne.tsv", "", basename(exp_file))
      output_dir <-
        file.path(output_base_dir, paste0(reg_name, "_", exp_file_base))
      create_directories(list(output_dir))
      execute_aracne(aracne_bin, exp_file, regulator_file, output_dir, threads,
                     seed)
    }
  }
}

#' Executes ARACNe3 on the provided expression matrix and regulator list.
#'
#' @param aracne_bin Path to the ARACNe3 binary.
#' @param exp_file Path to the expression matrix file.
#' @param regulators_file Path to the file containing regulator gene names.
#' @param output_dir Directory to store ARACNe output.
#' @param threads Number of threads to use for ARACNe computation.
#' @param seed Seed for random number generation in ARACNe.
execute_aracne <- function(aracne_bin, exp_file, regulators_file, output_dir,
                           threads = 1, seed = 123) {
  cmd <-
    sprintf("%s -e %s -r %s -o %s --threads %d --seed %d",
            aracne_bin, exp_file, regulators_file, output_dir, threads, seed)
  system(cmd)
}