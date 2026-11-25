source("classes/utils.r")

library(reshape2)
library(ggplot2)

#' MetacellGenerator Class for Generating Metacell Matrices
#'
#' This class provides methods to generate metacell matrices from a Seurat
#' object, which are the input for ARACNe network inference (see the `Runner`
#' class).
#'
#' A metacell is the sum of one cell's counts and the counts of its
#' `num_neighbors` nearest neighbours within the same cluster. One metacell is
#' built per cell, so neighbouring metacells share most of their cells. For
#' that reason the matrix used for ARACNe is a random subsample of `sub_size`
#' metacells per cluster, normalised to counts per million (CPM):
#'
#' - ARACNe3's significance test assumes independent samples; overlapping
#'   metacells violate that assumption and inflate false-positive edges.
#'   Subsampling keeps the overlap low for clusters much larger than
#'   `sub_size`.
#' - ARACNe3 requires input normalised for sequencing depth (its README, step
#'   1: "Normalize a gene expression profile for sequencing depth in each
#'   sample (CPM, TPM, etc.)").
#'
#' @field seurat_obj A Seurat object containing the single-cell data.
#' @field base_output_path A string specifying the base path where output
#'        files will be saved.
MetacellGenerator <- R6Class( # nolint
  "MetacellGenerator",
  public = list(
    seurat_obj = NULL,
    base_output_path = NULL,
    utils = Utils$new(),

    #' Initialize the MetacellGenerator Object
    #'
    #' @param seurat_obj       A Seurat object containing the single-cell data.
    #' @param base_output_path A string specifying the base path where output
    #'                         files will be saved.
    initialize = function(seurat_obj, base_output_path) {
      self$seurat_obj <- seurat_obj
      self$base_output_path <- base_output_path
    },

    #' Generate Metacell Matrices for ARACNe Analysis
    #'
    #' For every cluster with more than `size_thresh` cells, builds one
    #' metacell per cell and writes two versions to `base_output_path`:
    #'
    #' - `{out_name}_clust-{i}-metaCells_sub.{rds,tsv}`: the ARACNe input.
    #'   `sub_size` randomly chosen metacells (all of them if the cluster has
    #'   no more than `sub_size` cells), CPM-normalised. This is the matrix
    #'   returned for the cluster.
    #' - `{out_name}_clust-{i}-metaCells_all.{rds,tsv}`: every metacell, as
    #'   raw summed counts. Written for reference only; it is NOT an ARACNe
    #'   input (it is neither subsampled nor depth-normalised).
    #'
    #' `i` is the position of the cluster among the clusters that were kept
    #' (1, 2, ...), not its `seurat_clusters` label.
    #'
    #' @param out_name      Prefix for the saved metacell matrix files.
    #' @param size_thresh   Clusters with `size_thresh` cells or fewer are
    #'                      skipped.
    #' @param num_neighbors Number of nearest neighbours summed into each
    #'                      metacell (in addition to the cell itself).
    #' @param sub_size      Number of metacells kept per cluster for ARACNe,
    #'                      i.e. the number of samples ARACNe sees. Clusters
    #'                      with fewer cells keep all of their metacells.
    #' @param seed          Seed for the random subsample, so that reruns
    #'                      select the same metacells. The caller's random
    #'                      number generator state is restored afterwards. Use
    #'                      `NULL` to draw from the current state instead.
    #'
    #' @return              A named list of subsampled, CPM-normalised metacell
    #'                      matrices (genes x metacells), one per cluster, named
    #'                      `{out_name}_clust-{i}`. Pass it to `Runner$new()`.
    generate_metacell_matrices =
      function(out_name = "metacell", size_thresh = 50, num_neighbors = 5,
               sub_size = 200, seed = 42) {
        self$utils$check_whole_number(sub_size, "sub_size")
        if (!is.null(seed)) {
          self$utils$check_whole_number(seed, "seed", minimum = 0)
          restore_rng <- private$seed_rng(seed)
          on.exit(restore_rng(), add = TRUE)
        }

        counts_matrix <- as.matrix(self$seurat_obj[["SCT"]]@counts)
        clustering <- self$seurat_obj$seurat_clusters

        meta_mats <-
          private$create_metacell_matrices(
            counts_matrix, clustering, num_neighbors, sub_size,
            self$base_output_path, out_name, size_thresh
          )

        return(meta_mats)
      }
  ),

  #############################################################################
  #                           PRIVATE METHODS                                 #
  #############################################################################
  private = list(
    #' Seed the RNG and Return a Function That Undoes It
    #'
    #' Sets the seed and returns a function that restores the previous
    #' `.Random.seed` (or removes it if there was none). Callers register it
    #' with `on.exit()`, so seeding the metacell subsample has no side effect
    #' on the random numbers the user's session draws afterwards.
    #'
    #' @param seed The seed to set.
    #'
    #' @return     A function with no arguments that restores the RNG state.
    seed_rng = function(seed) {
      global <- globalenv()
      had_seed <- exists(".Random.seed", envir = global, inherits = FALSE)
      old_seed <- if (had_seed) get(".Random.seed", envir = global) else NULL
      set.seed(seed)
      function() {
        if (had_seed) {
          assign(".Random.seed", old_seed, envir = global)
        } else if (exists(".Random.seed", envir = global, inherits = FALSE)) {
          rm(".Random.seed", envir = global)
        }
      }
    },

    #' Create Metacell Matrices for ARACNe Analysis
    #'
    #' This function takes a gene expression matrix and its corresponding
    #' clustering, creates metacell matrices for each cluster, and saves these
    #' matrices. It is designed for use with ARACNe to facilitate analysis of
    #' gene regulatory networks.
    #'
    #' @param dat_mat       Matrix of raw gene expression (genes X samples).
    #' @param num_neighbors Number of neighbors to use for each metacell.
    #' @param clustering    Vector of cluster labels for each sample in
    #'                      `dat_mat`.
    #' @param sub_size      Number of metacells kept per cluster for ARACNe.
    #' @param out_dir       Directory where metacell matrices will be saved.
    #' @param out_name      Prefix for saved metacell matrix files.
    #' @param size_thresh   Minimum cluster size; clusters with this many cells
    #'                      or fewer will be ignored.
    #'
    #' @return              A named list of subsampled, CPM-normalised metacell
    #'                      matrices, one per cluster (see
    #'                      `generate_metacell_matrices()`).
    create_metacell_matrices =
      function(dat_mat, clustering, num_neighbors = 10, sub_size = 200,
               out_dir, out_name = "", size_thresh = 50) {
        if (is.null(dat_mat) || ncol(dat_mat) == 0) {
          stop("Input data matrix is empty or NULL.")
        }

        if (length(unique(clustering)) <= 1) {
          stop("Insufficient unique clusters for processing.")
        }

        # Generate cluster-specific matrices and filter out empty ones
        clust_mats <-
          private$generate_cluster_matrices(dat_mat, clustering, size_thresh)
        clust_mats <- private$filter_non_empty_matrices(clust_mats)

        # Initialize list to store metacell matrices
        meta_mats <- list()

        if (length(clust_mats) > 0) {
          for (i in seq_along(clust_mats)) {
            meta_mat <- private$process_cluster(clust_mats[[i]], num_neighbors,
                                                i, out_dir, out_name, sub_size)
            # The name identifies the cluster's files and, in Runner, its
            # ARACNe input, ARACNe output directory and regulon files.
            meta_mats[[paste0(out_name, "_clust-", i)]] <- meta_mat
          }
        } else {
          cat("No valid cluster matrices to process.\n")
        }

        return(meta_mats)
      },

    #' Generate Cluster-Specific Matrices
    #'
    #' This function divides a gene expression matrix into sub-matrices based
    #' on cluster labels, ensuring each sub-matrix contains only the data for a
    #' specific cluster. Clusters with a size below the specified threshold are
    #' ignored.
    #'
    #' @param dat_mat     Matrix of raw gene expression data, with genes as rows
    #'                    and samples as columns.
    #' @param clustering  A vector of cluster labels corresponding to each
    #'                    column in `dat_mat`.
    #' @param size_thresh Minimum size of clusters to be considered. Clusters
    #'                    smaller than this threshold will be ignored.
    #'
    #' @return            A list of matrices, each representing gene expression
    #'                    data for a specific cluster.
    generate_cluster_matrices = function(dat_mat, clustering, size_thresh) {
      private$cluster_matrices(dat_mat, clustering, size_thresh = size_thresh)
    },

    #' Generate and Optionally Save Cluster-Specific Matrices
    #'
    #' Splits the data matrix into cluster-specific matrices based on provided
    #' cluster labels. Can save the resulting matrices to files if a save path
    #' is specified.
    #'
    #' @param dat_mat     Data matrix to be split (features x samples).
    #' @param clust       Clustering labels for samples.
    #' @param save_path   Optional path for saving the resulting matrices; if
    #'                    not provided, matrices are returned in a list.
    #' @param save_pref   Optional prefix for file names when saving matrices.
    #' @param size_thresh Minimum number of samples required for a cluster to
    #'                    be processed; default is 100.
    #'
    #' @return            A list of matrices, one for each cluster, if
    #'                    `save_path` is not provided. Otherwise, nothing is
    #'                    explicitly returned.
    cluster_matrices = function(dat_mat, clust, save_path = NA,
                                save_pref = "", size_thresh = 100) {
      clust_table <- table(clust)
      clust_mats <- list()

      for (i in seq_along(clust_table)) {
        if (clust_table[i] > size_thresh) {
          clust_cells <- which(clust == names(clust_table)[i])

          # Check if clust_cells contains valid indices
          if (length(clust_cells) == 0 || max(clust_cells) > ncol(dat_mat)) {
            next
          }

          clust_mat <- dat_mat[, clust_cells, drop = FALSE]
          clust_mat <- clust_mat[rowSums(clust_mat) >= 1, , drop = FALSE]

          if (is.na(save_path)) {
            clust_mats[[names(clust_table)[i]]] <- clust_mat
          } else {
            private$save_cluster_matrix(clust_mat, save_path,
                                        save_pref, names(clust_table)[i])
          }
        }
      }

      if (is.na(save_path)) {
        return(clust_mats)
      }
    },

    #' Save a Cluster Matrix to an RDS File
    #'
    #' Saves a given cluster matrix to an RDS file, constructing the file name
    #' based on provided parameters.
    #'
    #' @param clust_mat  The cluster-specific matrix to save.
    #' @param save_path  Directory path where the file will be saved.
    #' @param save_pref  Prefix to be added to the file name.
    #' @param clust_name Name of the cluster, used in the file name.
    save_cluster_matrix = function(clust_mat, save_path, save_pref,
                                   clust_name) {
      file_path <-
        file.path(save_path, paste0(save_pref, "_", clust_name, ".rds"))
      saveRDS(clust_mat, file = file_path)
    },

    #' Filter Out Empty Cluster Matrices
    #'
    #' Removes any null entries from a list of matrices. This is typically used
    #' to exclude cluster-specific matrices that might have been deemed too
    #' small or otherwise invalid.
    #'
    #' @param clust_mats A list of matrices, where each matrix corresponds to a
    #'                   cluster's gene expression data.
    #'
    #' @return           A filtered list of matrices, with null entries
    #'                   removed.
    filter_non_empty_matrices = function(clust_mats) {
      clust_mats <- Filter(function(x) !is.null(x) && nrow(x) > 0, clust_mats)
      return(clust_mats)
    },

    #' Process Each Cluster to Generate Metacell Matrix
    #'
    #' Builds one metacell per cell of the cluster and saves two versions:
    #'
    #' - `{out_name}_clust-{cluster_idx}-metaCells_all`: every metacell, raw
    #'   summed counts. Reference only, not an ARACNe input.
    #' - `{out_name}_clust-{cluster_idx}-metaCells_sub`: `sub_size` randomly
    #'   chosen metacells (all of them if there are no more than `sub_size`),
    #'   CPM-normalised. This is the ARACNe input and the return value.
    #'
    #' @param mat           A matrix representing the gene expression data for
    #'                      a single cluster.
    #' @param num_neighbors The number of neighbors to consider for each
    #'                      metacell.
    #' @param cluster_idx   The index of the current cluster being processed.
    #' @param out_dir       The directory where output files will be saved.
    #' @param out_name      A prefix to be added to the names of the output
    #'                      files.
    #' @param sub_size      Number of metacells kept for ARACNe.
    #'
    #' @return              The subsampled, CPM-normalised metacell matrix,
    #'                      identical to the one saved as `..._sub`.
    process_cluster = function(mat, num_neighbors, cluster_idx, out_dir,
                               out_name, sub_size) {
      file_prefix <- paste0(out_name, "_clust-", cluster_idx, "-metaCells")

      # One metacell per cell; neighbouring metacells share most of their
      # cells. Duplicate gene names are dropped here, once, so the returned
      # matrix is exactly what is saved (ARACNe3 rejects duplicate genes).
      meta_mat <- private$meta_cells(mat, num_neighbors)
      meta_mat <- private$remove_duplicate_genes(meta_mat)

      private$save_meta_mat(meta_mat, out_dir, paste0(file_prefix, "_all"))

      # ARACNe input: a random subsample of metacells, which keeps their
      # overlap low, normalised for sequencing depth.
      if (sub_size < ncol(meta_mat)) {
        meta_mat <- meta_mat[, sample(colnames(meta_mat), sub_size),
                             drop = FALSE]
      }
      meta_mat <- private$cpmt_transform(meta_mat)

      private$save_meta_mat(meta_mat, out_dir, paste0(file_prefix, "_sub"))

      return(meta_mat)
    },

    #' Perform Counts Per Million (CPM) Normalization
    #'
    #' This function normalizes gene expression data using Counts Per Million
    #' (CPM) normalization, with an option for subsequent log2 transformation.
    #'
    #' @param dat_mat Matrix of gene expression data, with genes as rows and
    #'                samples as columns.
    #' @param l2      Logical indicating whether to apply log2 transformation
    #'                after CPM normalization.
    #'
    #' @return        A matrix with CPM-normalized (and optionally
    #'                log2-transformed) values.
    cpmt_transform = function(dat_mat, l2 = FALSE) {
      # Calculate CPM
      cpm_mat <- private$calculate_cpm(dat_mat)

      # Apply log2 transformation if specified
      if (l2) {
        cpm_mat <- private$log2_transform(cpm_mat)
      }

      return(cpm_mat)
    },

    #' Calculate Counts Per Million (CPM)
    #'
    #' Converts raw counts to CPM for normalization across samples.
    #'
    #' @param dat_mat Raw counts matrix with genes as rows and samples as
    #'                columns.
    #'
    #' @return        CPM-normalized matrix.
    calculate_cpm = function(dat_mat) {
      t(t(dat_mat) / colSums(dat_mat) * 1e6)
    },

    #' Apply Log2 Transformation
    #'
    #' Transforms CPM-normalized values using log2, adding 1 to avoid log of
    #' zero.
    #'
    #' @param cpm_mat Matrix of CPM-normalized gene expression values.
    #'
    #' @return        Log2-transformed matrix.
    log2_transform = function(cpm_mat) {
      log2(cpm_mat + 1)
    },

    #' Save Metacell Matrix as RDS and TSV
    #'
    #' Writes `{out_dir}/{file_prefix}.rds` and `{out_dir}/{file_prefix}.tsv`.
    #' Both files hold exactly the matrix passed in; no further subsampling
    #' happens here.
    #'
    #' @param meta_mat    The metacell matrix to be saved (genes x metacells).
    #' @param out_dir     The directory where the files will be saved.
    #' @param file_prefix The file name without extension.
    #'
    #' @return            None; the function's primary effect is to write
    #'                    files to disk.
    save_meta_mat = function(meta_mat, out_dir, file_prefix) {
      # Ensure the output directory exists
      if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
      }

      out_file <- file.path(out_dir, file_prefix)
      saveRDS(meta_mat, file = paste0(out_file, ".rds"))
      private$save_matrix_tsv(meta_mat, out_file)
    },

    #' Remove Duplicate Genes from the Matrix
    #'
    #' Keeps the first row for each gene name.
    #'
    #' @param dat_mat A matrix with genes as rows and samples as columns.
    #'
    #' @return        Matrix with duplicate genes removed.
    remove_duplicate_genes = function(dat_mat) {
      dat_mat[!duplicated(rownames(dat_mat)), , drop = FALSE]
    },

    #' Save Matrix as TSV File
    #'
    #' The first row holds "gene" followed by the metacell names; each further
    #' row holds a gene name followed by its values.
    #'
    #' @param dat_mat  A matrix with genes as rows and samples as columns.
    #' @param out_file Base path and name for the output TSV file.
    save_matrix_tsv = function(dat_mat, out_file) {
      formatted_matrix <- private$format_matrix_for_tsv(dat_mat)
      write.table(formatted_matrix, file = paste0(out_file, ".tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE,
                  col.names = FALSE)
    },

    #' Format Matrix for Saving as TSV
    #'
    #' @param dat_mat A matrix with genes as rows and samples as columns.
    #'
    #' @return        A matrix formatted for saving as TSV, including header
    #'                row.
    format_matrix_for_tsv = function(dat_mat) {
      sample_names <- colnames(dat_mat)
      gene_ids <- rownames(dat_mat)
      rbind(c("gene", sample_names), cbind(gene_ids, dat_mat))
    },

    #' Generate Meta Cell Matrix
    #'
    #' Creates a meta cell matrix by aggregating information from each cell's
    #' nearest neighbors. This can be useful for imputing missing values or
    #' enhancing signal in sparse datasets.
    #'
    #' @param dat_mat       A matrix of raw gene expression data
    #'                      (genes x samples).
    #' @param num_neighbors The number of nearest neighbors to consider for
    #'                      each cell.
    #' @param sub_size      Optional; if specified, subsets the resulting meta
    #'                      cell matrix to this number of cells.
    #'
    #' @return              A matrix representing meta cells, potentially
    #'                      subsetted.
    meta_cells = function(dat_mat, num_neighbors = 10, sub_size = NA) {
      if (num_neighbors >= ncol(dat_mat)) {
        stop("num_neighbors must be less than the number of columns in dat_mat")
      }
      dist_mat <- self$utils$compute_distance_matrix(dat_mat)
      knn_neighbors <- private$find_knn(dist_mat, num_neighbors)
      imp_mat <- private$impute_matrix(dat_mat, knn_neighbors)

      if (!is.na(sub_size) && sub_size > 0 && sub_size <= ncol(imp_mat)) {
        imp_mat <- private$subset_matrix(imp_mat, sub_size)
      }
      return(imp_mat)
    },

    #' Find K-Nearest Neighbors
    #'
    #' Identifies the k-nearest neighbors for each sample based on the
    #' distance matrix.
    #'
    #' @param dist_mat A distance matrix.
    #' @param k        The number of neighbors to identify.
    #'
    #' @return         A matrix indicating the indices of k-nearest neighbors
    #'                 for each sample.
    find_knn = function(dist_mat, k) {
      dist_mat <- as.matrix(dist_mat)

      knn_indices <- t(apply(dist_mat, 1, function(x) {
        actual_k <- min(k, length(x) - 1)
        order(x)[2:(actual_k + 1)]
      }))

      return(knn_indices)
    },

    #' Impute Matrix
    #'
    #' Creates an imputed matrix by aggregating the expression of each sample
    #' with its k-nearest neighbors.
    #'
    #' @param dat_mat       A matrix of gene expression data (genes x samples).
    #' @param knn_neighbors A matrix of k-nearest neighbor indices for each
    #'                      sample.
    #'
    #' @return              An imputed gene expression matrix.
    impute_matrix = function(dat_mat, knn_neighbors) {
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
    },

    #' Subset Matrix
    #'
    #' Subsets a matrix to a specified number of columns (samples), chosen
    #' randomly.
    #'
    #' @param mat      A matrix to be subsetted.
    #' @param sub_size The number of columns to retain in the subset.
    #'
    #' @return         A subsetted matrix.
    subset_matrix = function(mat, sub_size) {
      mat[, sample(ncol(mat), sub_size)]
    }
  )
)