source("utils.r")

library(reshape2)
library(ggplot2)

#' Generate Metacell Matrices for ARACNe Analysis
#'
#' @param seurat_obj   A Seurat object with integrated data.
#' @param out_dir      Output directory to save the metacell matrices.
#' @param out_name     Prefix for the saved metacell matrix files.
#' @param size_thresh  Minimum cluster size to consider.
#' @param num_neighbors Number of neighbors to use for each metacell.
#' @param sub_size     Target number of cells in each metacell.
#'
#' @return A list of metacell matrices, one per cluster.
generate_metacell_matrices <-
  function(seurat_obj, out_dir, out_name, size_thresh = 50, num_neighbors = 5,
           sub_size = 200) {
    counts_matrix <- as.matrix(seurat_obj[["SCT"]]@counts)
    clustering <- seurat_obj$seurat_clusters

    meta_mats <-
      create_metacell_matrices(counts_matrix, clustering, num_neighbors,
                               sub_size, out_dir, out_name, size_thresh)

    return(meta_mats)
  }

#' Create Metacell Matrices for ARACNe Analysis
#'
#' This function takes a gene expression matrix and its corresponding
#' clustering, creates metacell matrices for each cluster, and saves these
#' matrices. It is designed for use with ARACNe to facilitate analysis of gene
#' regulatory networks.
#'
#' @param dat_mat       Matrix of raw gene expression (genes X samples).
#' @param num_neighbors Number of neighbors to use for each metacell.
#' @param clustering    Vector of cluster labels for each sample in `dat_mat`.
#' @param sub_size      Target number of cells in each metacell for ARACNe
#'                      analysis.
#' @param out_dir       Directory where metacell matrices will be saved.
#' @param out_name      Prefix for saved metacell matrix files.
#' @param size_thresh   Minimum cluster size; clusters smaller than this will
#'                      be ignored.
#'
#' @return              A list of metacell matrices, one per cluster.
create_metacell_matrices <-
  function(dat_mat, clustering, num_neighbors = 10, sub_size = 200,
           out_dir, out_name = "", size_thresh = 50) {
    if (is.null(dat_mat) || ncol(dat_mat) == 0) {
      stop("Input data matrix is empty or NULL.")
    }

    if (length(unique(clustering)) <= 1) {
      stop("Insufficient unique clusters for processing.")
    }

    # Generate cluster-specific matrices and filter out empty ones
    clust_mats <- generate_cluster_matrices(dat_mat, clustering, size_thresh)
    clust_mats <- filter_non_empty_matrices(clust_mats)

    # Initialize list to store metacell matrices
    meta_mats <- list()

    if (length(clust_mats) > 0) {
      for (i in seq_along(clust_mats)) {
        meta_mat <- process_cluster(clust_mats[[i]], num_neighbors, i, out_dir,
                                    out_name, sub_size)
        meta_mats[[i]] <- meta_mat
      }
    } else {
      cat("No valid cluster matrices to process.\n")
    }

    return(meta_mats)
  }

#' Generate Cluster-Specific Matrices
#'
#' This function divides a gene expression matrix into sub-matrices based on
#' cluster labels, ensuring each sub-matrix contains only the data for a
#' specific cluster. Clusters with a size below the specified threshold are
#' ignored.
#'
#' @param dat_mat     Matrix of raw gene expression data, with genes as rows
#'                    and samples as columns.
#' @param clustering  A vector of cluster labels corresponding to each column
#'                    in `dat_mat`.
#' @param size_thresh Minimum size of clusters to be considered. Clusters
#'                    smaller than this threshold will be ignored.
#'
#' @return            A list of matrices, each representing gene expression
#'                    data for a specific cluster.
generate_cluster_matrices <- function(dat_mat, clustering, size_thresh) {
  cluster_matrices(dat_mat, clustering, size_thresh = size_thresh)
}

#' Generate and Optionally Save Cluster-Specific Matrices
#'
#' Splits the data matrix into cluster-specific matrices based on provided
#' cluster labels. Can save the resulting matrices to files if a save path is
#' specified.
#'
#' @param dat_mat     Data matrix to be split (features x samples).
#' @param clust       Clustering labels for samples.
#' @param save_path   Optional path for saving the resulting matrices; if not
#'                    provided, matrices are returned in a list.
#' @param save_pref   Optional prefix for file names when saving matrices.
#' @param size_thresh Minimum number of samples required for a cluster to be
#'                    processed; default is 100.
#'
#' @return            A list of matrices, one for each cluster, if `save_path`
#'                    is not provided. Otherwise, nothing is explicitly
#'                    returned.
cluster_matrices <- function(dat_mat, clust, save_path = NA,
                             save_pref = "", size_thresh = 100) {
  clust_table <- table(clust)
  clust_mats <- list()

  for (i in seq_along(clust_table)) {
    if (clust_table[i] > size_thresh) {
      clust_cells <- which(clust == names(clust_table)[i])

      # Check if clust_cells contains valid indices
      if (length(clust_cells) == 0 || max(clust_cells) > ncol(dat_mat)) {
        next  # Skip to the next iteration if no valid cells or out of bounds
      }

      clust_mat <- dat_mat[, clust_cells, drop = FALSE]
      clust_mat <- clust_mat[rowSums(clust_mat) >= 1, , drop = FALSE]

      if (is.na(save_path)) {
        clust_mats[[names(clust_table)[i]]] <- clust_mat
      } else {
        save_cluster_matrix(clust_mat, save_path,
                            save_pref, names(clust_table)[i])
      }
    }
  }

  if (is.na(save_path)) {
    return(clust_mats)
  }
}

#' Save a Cluster Matrix to an RDS File
#'
#' Saves a given cluster matrix to an RDS file, constructing the file name
#' based on provided parameters.
#'
#' @param clust_mat  The cluster-specific matrix to save.
#' @param save_path  Directory path where the file will be saved.
#' @param save_pref  Prefix to be added to the file name.
#' @param clust_name Name of the cluster, used in the file name.
save_cluster_matrix <- function(clust_mat, save_path, save_pref, clust_name) {
  file_path <- file.path(save_path, paste0(save_pref, "_", clust_name, ".rds"))
  saveRDS(clust_mat, file = file_path)
}

#' Filter Out Empty Cluster Matrices
#'
#' Removes any null entries from a list of matrices. This is typically used to
#' exclude cluster-specific matrices that might have been deemed too small or
#' otherwise invalid.
#'
#' @param clust_mats A list of matrices, where each matrix corresponds to a
#'                   cluster's gene expression data.
#'
#' @return           A filtered list of matrices, with null entries removed.
filter_non_empty_matrices <- function(clust_mats) {
  clust_mats <- Filter(function(x) !is.null(x) && nrow(x) > 0, clust_mats)
  return(clust_mats)
}

#' Process Each Cluster to Generate Metacell Matrix
#'
#' For a given cluster's gene expression matrix, this function generates a
#' metacell matrix by considering the specified number of neighbors. It saves
#' two versions of the metacell matrix: one with all cells and another with a
#' subset (if the original exceeds the `sub_size`). Both matrices are saved to
#' files.
#'
#' @param mat           A matrix representing the gene expression data for a
#'                      single cluster.
#' @param num_neighbors The number of neighbors to consider for each metacell.
#' @param cluster_idx   The index of the current cluster being processed.
#' @param out_dir       The directory where output files will be saved.
#' @param out_name      A prefix to be added to the names of the output files.
#' @param sub_size      The maximum number of cells to include in the subsetted
#'                      metacell matrix.
#'
#' @return              A metacell matrix for the cluster, potentially
#'                      subsetted and transformed.
process_cluster <- function(mat, num_neighbors, cluster_idx, out_dir,
                            out_name, sub_size) {
  # Generate metacell matrix
  meta_mat <- meta_cells(mat, num_neighbors)

  # Save the complete metacell matrix
  save_meta_mat(meta_mat, out_dir,
                paste0(out_name, "_clust-", cluster_idx, "-metaCells_all"),
                subset = FALSE)

  # Subset if necessary and apply CPM transformation
  if (sub_size < ncol(meta_mat)) {
    meta_mat <- meta_mat[, sample(colnames(meta_mat), sub_size)]
  }
  meta_mat <- cpmt_transform(meta_mat)

  # Save the subsetted and transformed metacell matrix
  save_meta_mat(meta_mat, out_dir,
                paste0(out_name, "_clust-", cluster_idx, "-metaCells"),
                subset = TRUE)

  return(meta_mat)
}

#' Perform Counts Per Million (CPM) Normalization
#'
#' This function normalizes gene expression data using Counts Per Million (CPM)
#' normalization, with an option for subsequent log2 transformation.
#'
#' @param dat_mat Matrix of gene expression data, with genes as rows and samples
#'                as columns.
#' @param l2      Logical indicating whether to apply log2 transformation after
#'                CPM normalization.
#'
#' @return        A matrix with CPM-normalized (and optionally
#'                log2-transformed) values.
cpmt_transform <- function(dat_mat, l2 = FALSE) {
  # Calculate CPM
  cpm_mat <- calculate_cpm(dat_mat)

  # Apply log2 transformation if specified
  if (l2) {
    cpm_mat <- log2_transform(cpm_mat)
  }

  return(cpm_mat)
}

#' Calculate Counts Per Million (CPM)
#'
#' Converts raw counts to CPM for normalization across samples.
#'
#' @param dat_mat Raw counts matrix with genes as rows and samples as columns.
#'
#' @return        CPM-normalized matrix.
calculate_cpm <- function(dat_mat) {
  t(t(dat_mat) / colSums(dat_mat) * 1e6)
}

#' Apply Log2 Transformation
#'
#' Transforms CPM-normalized values using log2, adding 1 to avoid log of zero.
#'
#' @param cpm_mat Matrix of CPM-normalized gene expression values.
#'
#' @return        Log2-transformed matrix.
log2_transform <- function(cpm_mat) {
  log2(cpm_mat + 1)
}

#' Save Metacell Matrix to File
#'
#' Saves a given metacell matrix to a file, constructing the file name from the
#' provided directory, file prefix, and an indicator of whether the matrix has
#' been subsetted.
#'
#' @param meta_mat    The metacell matrix to be saved.
#' @param out_dir     The directory where the file will be saved.
#' @param file_prefix The prefix to be used in constructing the file name.
#' @param subset      A boolean flag indicating whether the matrix is a
#'                    subsetted version.
#'
#' @return            None; the function's primary effect is to write a file
#'                    to disk.
save_meta_mat <- function(meta_mat, out_dir, file_prefix, subset) {
  # Ensure the output directory exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Construct the file name with appropriate suffix based on the subset flag
  file_suffix <- ifelse(subset, "_sub", "_all")
  file_name <- paste0(out_dir, "/", file_prefix, file_suffix, ".txt")

  # Call aracne_table to save the data appropriately
  aracne_table(meta_mat, file_name, subset)
}

#' Save Data Matrix for ARACNe Analysis
#'
#' Formats and saves a gene expression matrix for use as input to ARACNe,
#' optionally subsetting the matrix to a maximum number of samples for
#' efficiency.
#'
#' @param dat_mat  A matrix of data with genes as rows and samples as columns.
#' @param out_file Path and base name for the output file(s).
#' @param subset   Logical indicating whether to subset the matrix to 500
#'                 samples.
aracne_table <- function(dat_mat, out_file, subset = TRUE) {
  # Remove duplicate genes
  dat_mat <- remove_duplicate_genes(dat_mat)

  # Save the full or subsetted matrix as an RDS file
  save_matrix_rds(dat_mat, out_file, subset)

  # Prepare and save the matrix in TSV format
  save_matrix_tsv(dat_mat, out_file, subset)
}

#' Remove Duplicate Genes from the Matrix
#'
#' @param dat_mat A matrix with genes as rows and samples as columns.
#'
#' @return        Matrix with duplicate genes removed.
remove_duplicate_genes <- function(dat_mat) {
  dat_mat[!duplicated(rownames(dat_mat)), ]
}

#' Save Matrix as RDS File
#'
#' @param dat_mat  A matrix with genes as rows and samples as columns.
#' @param out_file Base path and name for the output RDS file.
#' @param subset   Logical indicating if the matrix should be subsetted.
save_matrix_rds <- function(dat_mat, out_file, subset) {
  if (subset) {
    dat_mat <- subset_matrix_samples(dat_mat, 500)
  }
  saveRDS(dat_mat, file = paste0(out_file, ".rds"))
}

#' Subset Matrix to a Specific Number of Samples
#'
#' @param dat_mat     A matrix with genes as rows and samples as columns.
#' @param max_samples The maximum number of samples to include in the subset.
#'
#' @return            A subsetted matrix.
subset_matrix_samples <- function(dat_mat, max_samples) {
  dat_mat[, sample(colnames(dat_mat), min(ncol(dat_mat), max_samples))]
}

#' Save Matrix as TSV File
#'
#' @param dat_mat  A matrix with genes as rows and samples as columns.
#' @param out_file Base path and name for the output TSV file.
#' @param subset   Logical indicating if the matrix should be subsetted before
#'                 saving.
save_matrix_tsv <- function(dat_mat, out_file, subset) {
  if (subset) {
    dat_mat <- subset_matrix_samples(dat_mat, 500)
  }
  formatted_matrix <- format_matrix_for_tsv(dat_mat)
  write.table(formatted_matrix, file = paste0(out_file, ".tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

#' Format Matrix for Saving as TSV
#'
#' @param dat_mat A matrix with genes as rows and samples as columns.
#'
#' @return        A matrix formatted for saving as TSV, including header row.
format_matrix_for_tsv <- function(dat_mat) {
  sample_names <- colnames(dat_mat)
  gene_ids <- rownames(dat_mat)
  rbind(c("gene", sample_names), cbind(gene_ids, dat_mat))
}

#' Generate Meta Cell Matrix
#'
#' Creates a meta cell matrix by aggregating information from each cell's
#' nearest neighbors. This can be useful for imputing missing values or
#' enhancing signal in sparse datasets.
#'
#' @param dat_mat       A matrix of raw gene expression data (genes x samples).
#' @param num_neighbors The number of nearest neighbors to consider for each
#'                      cell.
#' @param sub_size      Optional; if specified, subsets the resulting meta cell
#'                      matrix to this number of cells.
#'
#' @return              A matrix representing meta cells, potentially
#'                      subsetted.
meta_cells <- function(dat_mat, num_neighbors = 10, sub_size = NA) {
  if (num_neighbors >= ncol(dat_mat)) {
    stop("num_neighbors must be less than the number of columns in dat_mat")
  }
  # Compute distance matrix based on Pearson correlation
  dist_mat <- compute_distance_matrix(dat_mat)

  # Identify nearest neighbors for each sample
  knn_neighbors <- find_knn(dist_mat, num_neighbors)

  # Create imputed matrix based on nearest neighbors
  imp_mat <- impute_matrix(dat_mat, knn_neighbors)

  # Subset the imputed matrix if sub_size is specified and valid
  if (!is.na(sub_size) && sub_size > 0 && sub_size <= ncol(imp_mat)) {
    imp_mat <- subset_matrix(imp_mat, sub_size)
  }

  return(imp_mat)
}

#' Find K-Nearest Neighbors
#'
#' Identifies the k-nearest neighbors for each sample based on the distance
#' matrix.
#'
#' @param dist_mat A distance matrix.
#' @param k        The number of neighbors to identify.
#'
#' @return         A matrix indicating the indices of k-nearest neighbors for
#'                 each sample.
find_knn <- function(dist_mat, k) {
  # Ensure the distance matrix is in matrix format
  dist_mat <- as.matrix(dist_mat)

  # Apply over each row to find the indices of the k-nearest neighbors
  knn_indices <- t(apply(dist_mat, 1, function(x) {
    actual_k <- min(k, length(x) - 1)
    order(x)[2:(actual_k + 1)]
  }))

  return(knn_indices)
}

#' Impute Matrix
#'
#' Creates an imputed matrix by aggregating the expression of each sample with
#' its k-nearest neighbors.
#'
#' @param dat_mat       A matrix of gene expression data (genes x samples).
#' @param knn_neighbors A matrix of k-nearest neighbor indices for each sample.
#'
#' @return              An imputed gene expression matrix.
impute_matrix <- function(dat_mat, knn_neighbors) {
  # Initialize the imputed matrix
  imp_mat <- matrix(0, nrow = nrow(dat_mat), ncol = ncol(dat_mat))
  colnames(imp_mat) <- colnames(dat_mat)
  rownames(imp_mat) <- rownames(dat_mat)

  # Iterate over each sample to impute based on nearest neighbors
  for (i in seq_len(ncol(dat_mat))) {
    # Retrieve neighbor indices for the current sample
    neighbor_cols <- c(i, knn_neighbors[i, ])

    # Ensure all referenced indices are within bounds
    if (any(neighbor_cols > ncol(dat_mat))) {
      stop(paste("Out of bounds error at sample", i,
                 ": Neighbor indices",
                 toString(neighbor_cols[neighbor_cols > ncol(dat_mat)]),
                 "are greater than the number of columns", ncol(dat_mat)))
    }

    # Aggregate data from the original matrix using the neighbor indices
    imp_mat[, i] <-
      rowSums(dat_mat[, neighbor_cols, drop = FALSE], na.rm = TRUE)
  }

  return(imp_mat)
}

#' Subset Matrix
#'
#' Subsets a matrix to a specified number of columns (samples), chosen
#' randomly.
#'
#' @param mat      A matrix to be subsetted.
#' @param sub_size The number of columns to retain in the subset.
#'
#' @return         A subsetted matrix.
subset_matrix <- function(mat, sub_size) {
  mat[, sample(ncol(mat), sub_size)]
}

#' Plot Cluster Frequencies by Treatment
#'
#' @param seurat_obj       A Seurat object with integrated data.
#' @param plot_output_path Directory to save the plot.
#' @param col_names        Vector of column names for the combined data frame.
#' @param plot_title       Title of the plot.
#' @param binwidth         Width of bins in the dot plot.
#'
#' @return None. The function saves a plot to the specified directory.
plot_cluster_freq_by_treatment <-
  function(seurat_obj, plot_output_path,
           col_names = c("Early_p1", "Early_p2", "Early_p3", "Late_p1",
                         "Late_p2", "Late_p3"),
           plot_title = "Cluster Frequency by Early vs Late",
           binwidth = 0.01) {

    cluster_freq_table <-
      table(seurat_obj$patient, seurat_obj$seurat_clusters, seurat_obj$type)

    early_data <- calculate_cluster_frequencies(cluster_freq_table, 1)
    late_data <- calculate_cluster_frequencies(cluster_freq_table, 2)

    combined_data <-
      combine_cluster_frequencies(early_data, late_data, col_names)

    plot_data <- melt_cluster_frequencies(combined_data)

    plot_cluster_frequencies(plot_data, plot_title,
                             file.path(plot_output_path,
                                       "cluster_frequencies_by_treatment.png"),
                             binwidth)
  }

#' Calculate Cluster Frequencies
#'
#' @param cluster_freq_table Table of cluster frequencies.
#' @param treatment_index Index of the treatment type (1 for early, 2 for late).
#'
#' @return A data frame of normalized cluster frequencies.
calculate_cluster_frequencies <- function(cluster_freq_table, treatment_index) {
  freq_data <- as.data.frame.matrix(cluster_freq_table[, , treatment_index])
  freq_data <- freq_data[rowSums(freq_data) > 0, ]
  freq_data <- apply(freq_data, 1, function(x) {
    x / sum(x)
  })
  return(freq_data)
}

#' Combine Cluster Frequencies
#'
#' @param early_data Data frame of early treatment cluster frequencies.
#' @param late_data Data frame of late treatment cluster frequencies.
#' @param col_names Vector of column names for the combined data frame.
#'
#' @return A combined data frame of early and late treatment cluster
#'         frequencies.
combine_cluster_frequencies <- function(early_data, late_data, col_names) {
  combined <- cbind(early_data, late_data)
  colnames(combined) <- col_names
  return(combined)
}

#' Melt Cluster Frequencies
#'
#' @param combined_data Combined data frame of cluster frequencies.
#'
#' @return A melted data frame suitable for plotting.
melt_cluster_frequencies <- function(combined_data) {
  melted_data <- melt(combined_data)
  colnames(melted_data) <- c("cluster", "type", "frequency")
  melted_data$type <- unlist(
    lapply(strsplit(as.character(melted_data$type), "_"), function(x) {
      x[1]
    })
  )
  melted_data$type <- factor(melted_data$type, levels = c("Early", "Late"))
  melted_data$cluster <- as.factor(melted_data$cluster)
  return(melted_data)
}

#' Plot Cluster Frequencies
#'
#' @param plot_data Data frame of cluster frequencies to be plotted.
#' @param plot_title Title of the plot.
#' @param output_path File path to save the plot.
#' @param binwidth Width of bins in the dot plot.
#'
#' @return None. The function saves a plot to the specified file path.
plot_cluster_frequencies <-
  function(plot_data, plot_title, output_path, binwidth) {
    theme_update(plot.title = element_text(hjust = 0.5))

    p <- ggplot(plot_data, aes(x = cluster, y = frequency, fill = type)) + # nolint
      geom_dotplot(binaxis = "y", stackdir = "center",
                   position = position_dodge(), binwidth = binwidth) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")) +
      ggtitle(plot_title) +
      theme(plot.title = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 8),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 12),
            strip.text.x = element_text(size = 12, face = "bold"))

    ggsave(output_path, plot = p, width = 10, height = 8)
  }