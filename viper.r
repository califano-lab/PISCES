#' Process ARACNe Output Files to Generate Regulon Objects
#'
#' This function processes ARACNe output files to generate regulon objects
#' suitable for VIPER analysis. It loads ARACNe output files, prepares the
#' data, and generates regulon objects.
#'
#' @param aracne_output_base_dir The base directory containing ARACNe output
#'                               files.
#' @param exp_mat The expression matrix used to generate the ARACNe network
#'                (genes x samples).
#' @param output_base_path The base directory where the processed regulon
#'                         objects will be saved.
#'
#' @return A list of regulon objects generated from the ARACNe output files.
generate_regulon_objects <- function(aracne_output_base_dir, exp_mat,
                                     output_base_path) {
  aracne_output_files <- get_all_aracne_files(aracne_output_base_dir)

  if (length(aracne_output_files) == 0) {
    stop("No ARACNe output files found in the directory.")
  }

  regulon_list <- lapply(aracne_output_files, function(aracne_file) {
    aracne_data_for_viper <- prep_aracne_output_for_viper(aracne_file)
    prefix <- gsub("-metaCells$", "", basename(dirname(aracne_file)))
    regulon <- generate_regulon(aracne_data_for_viper, exp_mat,
                                output_base_path, prefix)
    return(regulon)
  })

  return(regulon_list)
}

#' Aggregate ARACNe Output Files from All Directories
#'
#' This helper function aggregates ARACNe output files from all directories
#' within the specified base directory, matching the provided file name
#' pattern.
#'
#' @param base_dir The base directory containing the ARACNe output files.
#' @param pattern A regular expression pattern to match ARACNe output files.
#'
#' @return A character vector of file paths to the ARACNe output files.
get_all_aracne_files <-
  function(base_dir, pattern = "consolidated-net_.*\\.tsv$") {
    all_files <- list.files(base_dir, pattern = pattern, full.names = TRUE,
                            recursive = TRUE)
    return(all_files)
  }

#' Load and Process ARACNe Output File
#'
#' This helper function loads an ARACNe output file, processes its content,
#' and prepares it for VIPER analysis.
#'
#' @param aracne_file The path to the ARACNe output file.
#'
#' @return A data frame containing the processed ARACNe data with columns
#'         "regulator", "target", and "mi" (mutual information).
prep_aracne_output_for_viper <- function(aracne_file) {
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
}

#' Generate Regulon Object from ARACNe Output
#'
#' This helper function generates a regulon object from ARACNe output data and
#' an expression matrix. The regulon object is then saved to the specified
#' output directory.
#'
#' @param aracne_data A data frame containing ARACNe output data.
#' @param exp_mat An expression matrix associated with the ARACNe network.
#' @param output_base_path The base directory where the regulon objects will
#'                         be saved.
#' @param file_prefix The prefix for the saved regulon files.
#'
#' @return A pruned regulon object suitable for VIPER analysis.
generate_regulon <- function(aracne_data, exp_mat, output_base_path,
                             file_prefix) {
  # Process ARACNe results for VIPER analysis
  reg_process(aracne_data, exp_mat, output_base_path, file_prefix)

  # Load the pruned regulon object
  pruned_regulon_file <-
    file.path(output_base_path, paste0(file_prefix, "_pruned.rds"))
  pruned_regulon <- readRDS(pruned_regulon_file)

  return(pruned_regulon)
}

#' Run VIPER on a List of Regulon Objects
#'
#' This function runs VIPER analysis on a list of regulon objects using an
#' expression matrix. It processes each regulon object and returns the VIPER
#' scores.
#'
#' @param exp_mat An expression matrix used for VIPER analysis.
#' @param regulon_list A list of regulon objects generated from ARACNe output.
#'
#' @return A list of VIPER results for each regulon object.
run_viper <- function(exp_mat, regulon_list) {
  viper_results <- lapply(regulon_list, function(regulon) {
    viper_scores <- execute_viper(exp_mat, regulon)
    return(viper_scores)
  })

  # Check contents of viper_results
  if (length(viper_results) == 0 || any(sapply(viper_results, is.null))) {
    stop("VIPER results are empty or not properly formed.")
  }

  return(viper_results)
}

#' Process ARACNe Results for VIPER Analysis
#'
#' Converts ARACNe output into a regulon object suitable for VIPER analysis,
#' including an optional pruning step to refine the regulon.
#'
#' @param a_file   Path to the ARACNe final network file in TSV format.
#' @param exp_mat  Expression matrix used to generate the ARACNe network
#'                 (genes x samples).
#' @param out_dir  Directory where the processed regulon objects will be saved.
#' @param out_name Prefix for the saved regulon files.
reg_process <- function(a_file, exp_mat, out_dir, out_name = "") {
  require(viper)

  # Convert ARACNe output to regulon object
  processed_reg <- convert_to_regulon(a_file, exp_mat)

  # Save the unpruned regulon object
  save_regulon(processed_reg, out_dir, paste0(out_name, "_unpruned.rds"))

  # Prune the regulon to refine it
  pruned_reg <- prune_regulon(processed_reg)

  # Save the pruned regulon object
  save_regulon(pruned_reg, out_dir, paste0(out_name, "_pruned.rds"))
}

#' Convert ARACNe Output to Regulon Object
#'
#' @param a_file  Path to the ARACNe network file.
#' @param exp_mat Expression matrix associated with the ARACNe network.
#'
#' @return        A regulon object suitable for VIPER analysis.
convert_to_regulon <- function(aracne_data, exp_mat) {
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
}

#' Save Regulon Object to File
#'
#' @param regulon   The regulon object to be saved.
#' @param out_dir   The directory for saving the regulon file.
#' @param file_name The name of the file to save the regulon object in.
save_regulon <- function(regulon, out_dir, file_name) {
  saveRDS(regulon, file = file.path(out_dir, file_name))
}

#' Prune Regulon Object
#'
#' Applies pruning to a regulon object to refine its content.
#'
#' @param regulon The regulon object to be pruned.
#'
#' @return        A pruned regulon object.
prune_regulon <- function(regulon) {
  viper::pruneRegulon(regulon, 50, adaptive = FALSE, eliminate = TRUE)
}

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
execute_viper <- function(exp_mat, regulon) {
  viper_scores <- tryCatch({
    viper(exp_mat, regulon)
  }, error = function(e) {
    cat("Error in VIPER analysis:", e$message, "\n")
    NULL
  })

  return(viper_scores)
}

#' Save VIPER Results to File
#'
#' This function saves the VIPER analysis results to a specified file in
#' RDS format.
#'
#' @param viper_results A list of VIPER results to be saved.
#' @param output_path The file path where the VIPER results will be saved.
#'
#' @return None.
save_viper_results <- function(viper_results, output_path) {
  saveRDS(viper_results, file = output_path)
}