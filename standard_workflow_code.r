library(Atools)
library(MASS)
library(reticulate)
library(SingleR)
library(Seurat)
library(cowplot)
library(dplyr)
library(cluster)
library(umap)
library(reshape)
library(pheatmap)
library(viper)
library("org.Hs.eg.db")
library(clustree)
library(factoextra)
library(MAST)
library(Hmisc)
library(ggplot2)
library(scales)
library(flowCore)
library(ggcyto)
library(infercnv)
library(ggrepel)
library(plyr)
library(PISCES)

#' Plot a heatmap of custom gene list grouped by cluster
#'
#' This function plots a heatmap for a subset of genes, potentially grouped by
#' cluster. It allows for customization of the color palette, scaling, and
#' selection of top genes per cluster.
#'
#' @param dat                     Matrix with genes as rows and samples as
#'                                columns.
#' @param clust                   Vector of cluster labels.
#' @param genes                   Vector of genes to include in the heatmap.
#' @param genes_by_cluster        Whether to group genes by cluster identity.
#' @param n_top_genes_per_cluster Number of top genes per cluster to plot.
#' @param color_palette           Custom color palette for clusters; uses
#'                                hue_pal by default.
#' @param scaled                  Whether the data is already scaled; applies
#'                                row-wise z-score scaling if FALSE.
#'
#' @return                        Heatmap plot.
gene_heatmap_plot <- function(dat, clust, genes, genes_by_cluster = TRUE,
                              n_top_genes_per_cluster = 5, color_palette = NA,
                              scaled = FALSE) {
  identities <- levels(clust)
  my_color_palette <- generate_color_palette(identities, color_palette)

  i <- sample(seq_len(ncol(dat)), min(10000, ncol(dat)), replace = FALSE)
  x <- dat[genes, i]

  df <- data.frame(cluster = clust[i])
  rownames(df) <- colnames(x)

  o <- order(df$cluster)
  x <- x[, o]
  df <- df[o, , drop = FALSE]

  if (!scaled) {
    t <- apply(x, 1, calculate_z_score)
  } else {
    t <- x
  }

  mat_breaks <- generate_mat_breaks(t)
  annotations <- generate_annotations(df, my_color_palette, genes_by_cluster,
                                      n_top_genes_per_cluster)

  pheatmap_args <- list(x, cluster_rows = FALSE, show_rownames = TRUE,
                        cluster_cols = FALSE, annotation_col = df,
                        breaks = mat_breaks,
                        color = colorRampPalette(c("blue",
                                                   "white",
                                                   "red"))
                        (length(mat_breaks)),
                        fontsize_row = ifelse(genes_by_cluster, 10, 8),
                        show_colnames = FALSE,
                        annotation_colors = annotations$anno_colors)

  if (!is.null(annotations$anno_row)) {
    pheatmap_args$annotation_row <- annotations$anno_row
    pheatmap_args$gaps_row <-
      (2:length(unique(clust)) - 1) * n_top_genes_per_cluster
  }

  do.call(pheatmap, pheatmap_args)
}

#' Generate a color palette
#'
#' Generates a color palette based on the provided identities. Defaults to
#' hue_pal if no custom palette is provided.
#'
#' @param identities    Vector of unique identities for which colors are
#'                      needed.
#' @param color_palette Optional custom color palette.
#'
#' @return              A color palette vector.
generate_color_palette <- function(identities, color_palette) {
  if (is.na(color_palette)) {
    return(hue_pal()(length(identities)))
  } else {
    return(color_palette)
  }
}

#' Calculate row-wise z-score
#'
#' Applies z-score normalization across rows of a matrix.
#'
#' @param x A numeric matrix.
#' @return  Matrix with row-wise z-scores.
calculate_z_score <- function(x) {
  return((x - mean(x)) / sd(x))
}

#' Generate matrix breaks based on quantiles
#'
#' Determines breaks for the heatmap color scale based on quantiles, excluding
#' extreme values.
#'
#' @param t Numeric matrix for which to determine breaks.
#'
#' @return  Vector of breaks for the heatmap color scale.
generate_mat_breaks <- function(t) {
  quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n + 1))
    breaks[!duplicated(breaks)]
  }

  lower_breaks <- quantile_breaks(t[t < 0], n = 10)
  upper_breaks <- quantile_breaks(t[t > 0], n = 10)
  c(lower_breaks, 0, upper_breaks)[-c(1, length(lower_breaks),
                                      length(lower_breaks) + 2,
                                      length(lower_breaks) +
                                        length(upper_breaks) + 1)]
}


#' Generate annotations for heatmap
#'
#' Creates a list containing colors for cluster annotations and an optional
#' data frame for row annotations if genes are grouped by cluster. The colors
#' are matched to clusters, and if genes are grouped by cluster, each gene
#' group is annotated with its corresponding cluster.
#'
#' @param df                      Data frame containing cluster information for
#'                                columns in the heatmap.
#' @param my_color_palette        Vector of colors used for cluster
#'                                annotations.
#' @param genes_by_cluster        Boolean indicating whether genes should be
#'                                grouped by their cluster.
#' @param n_top_genes_per_cluster Number of top genes per cluster to include if
#'                                genes are grouped by cluster.
#'
#' @return                        A list containing 'anno_colors' for column
#'                                annotations and 'anno_row' for row
#'                                annotations (if applicable).
generate_annotations <- function(df, my_color_palette, genes_by_cluster,
                                 n_top_genes_per_cluster) {
  anno_colors <- list(cluster = my_color_palette)
  names(anno_colors$cluster) <- levels(df$cluster)

  if (genes_by_cluster) {
    anno_colors$group <- anno_colors$cluster
    anno_row <- data.frame(group = rep(levels(df$cluster),
                                       each = n_top_genes_per_cluster))
    return(list(anno_colors = anno_colors, anno_row = anno_row))
  } else {
    return(list(anno_colors = anno_colors, anno_row = NULL))
  }
}


#' Compute silhouette scores for multiple Louvain clusterings
#'
#' This function evaluates multiple Louvain clusterings of a single-cell matrix
#' by computing the silhouette scores for a subsample of the data. It iterates
#' over 100 alternative resolution values, each time sub-sampling 1000 cells
#' (or fewer if fewer are available) and calculates the mean and standard
#' deviation of the silhouette scores for each clustering resolution.
#'
#' @param mat   Matrix with rows as principal component vectors and columns as
#'              samples.
#' @param clust Matrix with rows as samples and each column as a clustering
#'              vector for a given resolution.
#'
#' @return      List containing the means and standard deviations of silhouette
#'              scores for each clustering resolution.
sil_subsample <- function(mat, clust) {
  num_resolutions <- ncol(clust)
  num_subsamples <- 100

  # Initialize output dataframe
  silhouette_scores <- initialize_silhouette_scores(num_subsamples,
                                                    num_resolutions)

  for (resolution_index in 1:num_resolutions) {
    for (subsample_index in 1:num_subsamples) {
      sampled_indices <- sample_cells(mat, 1000)
      distance_matrix <- compute_distance_matrix(mat[, sampled_indices])
      silhouette_scores[subsample_index, resolution_index] <-
        compute_silhouette_width(clust[sampled_indices, resolution_index],
                                 distance_matrix)
    }
  }

  list(means = colMeans(silhouette_scores, na.rm = TRUE),
       sd = apply(silhouette_scores, 2, sd, na.rm = TRUE))
}

#' Initialize matrix for storing silhouette scores
#'
#' @param num_subsamples  Number of subsamples to compute.
#' @param num_resolutions Number of resolution values/clustering vectors.
#'
#' @return                Initialized matrix for storing silhouette scores.
initialize_silhouette_scores <- function(num_subsamples, num_resolutions) {
  matrix(rep(NA, num_subsamples * num_resolutions), nrow = num_subsamples)
}

#' Randomly sample cells from the matrix
#'
#' @param mat       Data matrix.
#' @param num_cells Number of cells to sample.
#'
#' @return          Indices of sampled cells.
sample_cells <- function(mat, num_cells) {
  sample(seq_len(ncol(mat)), min(num_cells, ncol(mat)))
}

#' Compute distance matrix using Pearson correlation
#'
#' @param mat Subsampled data matrix.
#'
#' @return Distance matrix.
compute_distance_matrix <- function(mat) {
  as.dist(1 - cor(mat, method = "pearson"))
}

#' Compute the mean silhouette width for a clustering
#'
#' @param clustering      Clustering vector for the subsampled cells.
#' @param distance_matrix Distance matrix for the subsampled cells.
#'
#' @return                Mean silhouette width for the clustering.
compute_silhouette_width <- function(clustering, distance_matrix) {
  if (length(unique(clustering)) <= 1) return(0)

  silhouette_scores <- silhouette(as.numeric(clustering), distance_matrix)
  mean(silhouette_scores[, "sil_width"])
}


#' Create a large SingleR object with improved control over tuning
#'
#' This function creates a SingleR object for large datasets by processing in
#' chunks, addressing the issue in the original SingleR where fine tuning was
#' always enabled.
#'
#' @param counts                Expression count matrix with genes in rows and
#'                              cells in columns.
#' @param annot                 Optional annotation vector for cells.
#' @param project.name          Name of the project for file naming.
#' @param xy                    Coordinates for cells, typically used for
#'                              visualization.
#' @param clusters              Cluster assignments for cells.
#' @param N                     Number of cells to process in each chunk.
#' @param min.genes             Minimum number of genes for filtering cells.
#' @param technology            Single-cell technology used (e.g., "10X").
#' @param species               Species name (e.g., "Human").
#' @param citation              Citation or reference for the dataset.
#' @param ref.list              Reference list for SingleR classification.
#' @param normalize.gene.length Logical indicating whether to normalize by gene
#'                              length.
#' @param variable.genes        Method for selecting variable genes, default is
#'                              "de" (differentially expressed).
#' @param fine.tune             Logical indicating whether to fine-tune SingleR
#'                              results.
#' @param reduce.file.size      Logical indicating whether to reduce file size
#'                              by omitting some SingleR slots.
#' @param do.signatures         Logical indicating whether to compute signature
#'                              scores.
#' @param do.main.types         Logical indicating whether to classify main
#'                              cell types.
#' @param temp.dir              Directory to store temporary SingleR objects.
#' @param numCores              Number of cores to use for parallel processing.
#'
#' @return                      A combined SingleR object created from chunks.
create_big_single_r_object_v2 <- function(counts, annot = NULL, project_name,
                                          xy, clusters, n = 10000,
                                          min_genes = 200, technology = "10X",
                                          species = "Human", citation = "",
                                          ref_list = list(),
                                          normalize_gene_length = FALSE,
                                          variable_genes = "de",
                                          fine_tune = TRUE,
                                          reduce_file_size = TRUE,
                                          do_signatures = FALSE,
                                          do_main_types = TRUE,
                                          temp_dir = getwd(), num_cores = 1) {
  setup_temp_dir(temp_dir, project_name)
  cell_indices <- split_into_chunks(ncol(counts), n)

  for (A in cell_indices) {
    process_chunk(A, counts, annot, project_name, min_genes, technology,
                  species, citation, do_signatures, num_cores, fine_tune,
                  temp_dir)
  }

  singler_objects <- load_singler_objects(temp_dir, project_name)
  combine_singler_objects(singler_objects, colnames(counts), clusters, xy)
}

#' Setup Temporary Directory for SingleR Objects
#'
#' Creates a temporary directory within the specified path to store SingleR
#' object chunks.
#'
#' @param temp_dir     The base directory to create a temporary directory in.
#' @param project_name The name of the project, used in naming the temporary
#'                     directory.
setup_temp_dir <- function(temp_dir, project_name) {
  dir.create(file.path(temp_dir, "singler_temp"), showWarnings = FALSE)
}

#' Split Total Cells into Chunks
#'
#' Divides the total number of cells into smaller chunks for processing.
#'
#' @param total_cells Total number of cells in the dataset.
#' @param chunk_size  Desired number of cells in each chunk.
#'
#' @return            A list where each element contains cell indices for a
#'                    chunk.
split_into_chunks <- function(total_cells, chunk_size) {
  split(seq_len(total_cells), ceiling(seq_len(total_cells) / chunk_size))
}

#' Process a Chunk of Cells
#'
#' Processes a subset of cells to create a SingleR object for that chunk.
#'
#' @param cell_indices Indices of cells in the chunk.
#' @param counts Expression counts matrix.
#' @param annot Cell annotations.
#' @param project_name Name of the project.
#' @param min_genes Minimum number of genes for inclusion.
#' @param technology Single-cell sequencing technology used.
#' @param species Species of the samples.
#' @param citation Citation for the dataset.
#' @param do_signatures Whether to compute signature scores.
#' @param num_cores Number of cores to use for computation.
#' @param fine_tune Whether to fine-tune the SingleR results.
#' @param temp_dir Temporary directory for storing SingleR objects.
process_chunk <- function(cell_indices, counts, annot, project_name, min_genes,
                          technology, species, citation, do_signatures,
                          num_cores, fine_tune, temp_dir) {
  singler <- SingleR::CreateSinglerObject(counts[, cell_indices],
                                          annot = annot[cell_indices],
                                          project_name = project_name,
                                          min_genes = min_genes,
                                          technology = technology,
                                          species = species,
                                          citation = citation,
                                          do_signatures = do_signatures,
                                          clusters = NULL,
                                          num_cores = num_cores,
                                          fine_tune = fine_tune)
  save(singler, file = file.path(temp_dir, "singler_temp",
                                 paste0(project_name, ".",
                                        cell_indices[1], ".RData")))
}

#' Load SingleR Objects from Files
#'
#' Loads SingleR objects saved in temporary files into a list.
#'
#' @param temp_dir     Directory where SingleR temporary files are stored.
#' @param project_name Name of the project, used to identify relevant files.
#'
#' @return             A list of SingleR objects.
load_singler_objects <- function(temp_dir, project_name) {
  singler_files <- list.files(file.path(temp_dir, "singler_temp"),
                              pattern = "RData", full.names = TRUE)

  lapply(singler_files, function(f) {
    # Load the .RData file
    loaded_names <- load(f)

    # Assume the object of interest is the last one loaded
    last_loaded_name <- tail(loaded_names, n = 1)

    # Use get() to retrieve the last loaded object by name
    get(last_loaded_name, envir = .GlobalEnv)
  })
}


#' Combine Singler objects into one
#'
#' @param singler_objects A list of Singler objects to combine.
#' @param cell_order      The order of cells to be maintained in the combined
#'                        object.
#' @param clusters        Cluster assignments for the cells.
#' @param xy              Coordinates for the cells, typically used for
#'                        visualization.
#'
#' @return                A combined Singler object.
combine_singler_objects <- function(singler_objects, cell_order,
                                    clusters, xy) {
  SingleR::SingleR.Combine(singler_objects, order = cell_order,
                           clusters = clusters, xy = xy)
}


#' Identify and Merge Unique Master Regulators for Each Cell
#'
#' This function computes the top master regulators (MRs) based on protein
#' activity for each cell in the input dataset. It then merges these lists and
#' returns a unique set of MRs across all cells.
#'
#' @param dat_mat A matrix with proteins as rows and samples (cells) as
#'                columns.
#' @param num_mrs Number of top MRs to identify in each cell.
#'
#' @return        A vector of unique master regulators identified across all
#'                cells.
cbcmrs <- function(dat_mat, num_mrs = 25) {
  # Identify MRs for each cell and merge into a unique list
  unique_mrs <- identify_and_merge_mrs(dat_mat, num_mrs)
  return(unique_mrs)
}

#' Identify and Merge MRs from Protein Activity Matrix
#'
#' This helper function applies over the columns of the protein activity matrix
#' to identify top MRs for each cell and then merges these lists, ensuring
#' uniqueness.
#'
#' @param dat_mat A matrix with proteins as rows and samples (cells) as
#'                columns.
#' @param num_mrs Number of top MRs to identify in each cell.
#'
#' @return        A vector of unique master regulators identified across all
#'                cells.
identify_and_merge_mrs <- function(dat_mat, num_mrs) {
  cbc_mrs <- apply(dat_mat, 2, function(x) {
    names(sort(x, decreasing = TRUE))[1:num_mrs]
  })
  unique(unlist(cbc_mrs))
}

#' Generate Metacell Matrices for ARACNe Analysis
#'
#' This function takes a gene expression matrix and its corresponding
#' clustering, generates metacell matrices for each cluster, and saves these
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
make_cmfa <- function(dat_mat, clustering, num_neighbors = 10, sub_size = 200,
                      out_dir, out_name = "", size_thresh = 100) {
  # Generate cluster-specific matrices and filter out empty ones
  clust_mats <- generate_cluster_matrices(dat_mat, clustering, size_thresh)
  clust_mats <- filter_non_empty_matrices(clust_mats)

  # Initialize list to store metacell matrices
  meta_mats <- list()

  for (i in seq_along(clust_mats)) {
    # Generate and save metacell matrix for each cluster
    meta_mat <- process_cluster(clust_mats[[i]], num_neighbors, i, out_dir,
                                out_name, sub_size)
    meta_mats[[i]] <- meta_mat
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
  clust_mats[which(!unlist(lapply(clust_mats, is.null)))]
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
  meta_mat <- CPMTransform(meta_mat)

  # Save the subsetted and transformed metacell matrix
  save_meta_mat(meta_mat, out_dir,
                paste0(out_name, "_clust-", cluster_idx, "-metaCells"),
                subset = TRUE)

  return(meta_mat)
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
  file_name <- paste0(out_dir, "/", file_prefix,
                      ifelse(subset, "", "_sub"), ".txt")
  ARACNeTable(meta_mat, file_name, subset)
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

#' Compute Distance Matrix
#'
#' Calculates a distance matrix using Pearson correlation.
#'
#' @param dat_mat A matrix of gene expression data (genes x samples).
#'
#' @return        A distance matrix.
compute_distance_matrix <- function(dat_mat) {
  as.dist(1 - cor(dat_mat, method = "pearson"))
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
  apply(as.matrix(dist_mat), 1, function(x) order(x)[2:(k + 1)])
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
  imp_mat <- matrix(0, nrow = nrow(dat_mat), ncol = ncol(dat_mat))
  colnames(imp_mat) <- colnames(dat_mat)
  rownames(imp_mat) <- rownames(dat_mat)

  for (i in seq_len(ncol(dat_mat))) {
    neighbor_cols <- c(i, knn_neighbors[i, ])
    imp_mat[, i] <- rowSums(dat_mat[, neighbor_cols, drop = FALSE])
  }

  imp_mat
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
      clust_mat <- dat_mat[, clust_cells]
      clust_mat <- clust_mat[rowSums(clust_mat) >= 1, ]

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

#' Performs a CPM normalization on the given data.
#'
#' @param dat.mat Matrix of gene expression data (genes X samples).
#' @param l2 Optional log2 normalization switch. Default of False.
#' @return Returns CPM normalized matrix
CPMTransform <- function(dat.mat, l2 = FALSE) {
  cpm.mat <- t(t(dat.mat) / (colSums(dat.mat) / 1e6))
  if (l2) {
    cpm.mat <- log2(cpm.mat + 1)
  }
  return(cpm.mat)
}

#' Saves a matrix in a format for input to ARACNe
#'
#' @param dat.mat Matrix of data (genes X samples).
#' @param out.file Output file where matrix will be saved.
#' @param subset Switch for subsetting the matrix to 500 samples. Default TRUE.
ARACNeTable <- function(dat.mat, out.file, subset = TRUE) {
  dat.mat <- dat.mat[!duplicated(rownames(dat.mat)), ]
  saveRDS(dat.mat, file = paste(out.file, ".rds", sep = ""))
  if (subset) {
    dat.mat <- dat.mat[, sample(colnames(dat.mat), min(ncol(dat.mat), 500))]
  }
  sample.names <- colnames(dat.mat)
  gene.ids <- rownames(dat.mat)
  m <- dat.mat
  mm <- rbind(c("gene", sample.names), cbind(gene.ids, m))
  write.table(
    x = mm, file = paste(out.file, ".tsv", sep = ""),
    sep = "\t", quote = F, row.names = F, col.names = F
  )
}

#' Processes ARACNe results into a regulon object compatible with VIPER.
#'
#' @param a.file ARACNe final network .tsv.
#' @param exp.mat Matrix of expression from which the network was generated (genes X samples).
#' @param out.dir Output directory for networks to be saved to.
#' @param out.name Optional argument for prefix of the file name.
RegProcess <- function(a.file, exp.mat, out.dir, out.name = ".") {
  require(viper)
  processed.reg <- aracne2regulon(afile = a.file, eset = exp.mat, format = "3col")
  saveRDS(processed.reg, file = paste(out.dir, out.name, "unPruned.rds", sep = ""))
  pruned.reg <- pruneRegulon(processed.reg, 50, adaptive = FALSE, eliminate = TRUE)
  saveRDS(pruned.reg, file = paste(out.dir, out.name, "pruned.rds", sep = ""))
}

fishersMethod <- function(x) {
  return(pchisq(-2 * sum(log(x)), df = 2 * length(x), lower = F))
} # integrate p-values
stouffersMethod <- function(x, weights = F) {
  if (weights) {
    return(sum(x * weights) / sqrt(sum(weights * weights)))
  } else {
    return(sum(x) / sqrt(length(x)))
  }
} # integrate z-scores



##### load raw data
neoredp1 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/victor_nygc_samples/neo-RED-P-C001-HIMC-CD45pos-GEX/outs/filtered_feature_bc_matrix/")
neoredp1 <- CreateSeuratObject(counts = neoredp1, min.features = 200, min.cells = 50)
neoredp1$patient <- "Patient1"
neoredp1$treatment <- "ADT+aCTLA4"
neoredp1$tissue <- "CD45+"
neoredp2 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/victor_nygc_samples/neo-RED-P-C001-HIMC-Total-GEX/outs/filtered_feature_bc_matrix/")
neoredp2 <- CreateSeuratObject(counts = neoredp2, min.features = 200, min.cells = 50)
neoredp2$patient <- "Patient1"
neoredp2$treatment <- "ADT+aCTLA4"
neoredp2$tissue <- "Total"
neoredp3 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/victor_nygc_samples/neo-RED-P-C002-Miltenyi-CD45pos-GEX/outs/filtered_feature_bc_matrix/")
neoredp3 <- CreateSeuratObject(counts = neoredp3, min.features = 200, min.cells = 50)
neoredp3$patient <- "Patient2"
neoredp3$treatment <- "ADT+aCTLA4"
neoredp3$tissue <- "CD45+"
neoredp4 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/victor_nygc_samples/neo-RED-P-C002-Miltenyi-Total-GEX/outs/filtered_feature_bc_matrix/")
neoredp4 <- CreateSeuratObject(counts = neoredp4, min.features = 200, min.cells = 50)
neoredp4$patient <- "Patient2"
neoredp4$treatment <- "ADT+aCTLA4"
neoredp4$tissue <- "Total"
# neoredp5<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/victor_nygc_samples/MC08-neo-RED-P-P003-U-no-treatment-total-GEX/filtered_feature_bc_matrix/")
# neoredp5<-CreateSeuratObject(counts = neoredp5,min.features = 200,min.cells=50)
# neoredp5$patient="Patient3"
# neoredp5$treatment="Untreated"
# neoredp5$tissue="Total"
# neoredp6<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/victor_nygc_samples/MC09-neo-RED-P-P003-U-no-treatment-CD45pos-GEX/filtered_feature_bc_matrix/")
# neoredp6<-CreateSeuratObject(counts = neoredp6,min.features = 200,min.cells=50)
# neoredp6$patient="Patient3"
# neoredp6$treatment="Untreated"
# neoredp6$tissue="CD45+"
neoredp7 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/MC11-neo-RED-P-P004-C-CTLA4-deg-CD45-GEX/filtered_feature_bc_matrix")
neoredp7 <- CreateSeuratObject(counts = neoredp7, min.features = 200, min.cells = 50)
neoredp7$patient <- "Patient4"
neoredp7$treatment <- "ADT+aCTLA4"
neoredp7$tissue <- "CD45+"
neoredp8 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/MC10-neo-RED-P-P004-C-CTLA4-deg-total-GEX/filtered_feature_bc_matrix")
neoredp8 <- CreateSeuratObject(counts = neoredp8, min.features = 200, min.cells = 50)
neoredp8$patient <- "Patient4"
neoredp8$treatment <- "ADT+aCTLA4"
neoredp8$tissue <- "Total"
# neoredp9<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/MC13-neo-RED-P-P005-U-no-treatment-CD45-GEX/filtered_feature_bc_matrix")
# neoredp9<-CreateSeuratObject(counts = neoredp9,min.features = 200,min.cells=50)
# neoredp9$patient="Patient5"
# neoredp9$treatment="Untreated"
# neoredp9$tissue="CD45+"
# neoredp10<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/MC12-neo-RED-P-P005-U-no-treatment-total-GEX/filtered_feature_bc_matrix")
# neoredp10<-CreateSeuratObject(counts = neoredp10,min.features = 200,min.cells=50)
# neoredp10$patient="Patient5"
# neoredp10$treatment="Untreated"
# neoredp10$tissue="Total"
neoredp11 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/MC15_neo-RED-P_P005-C_CTLA4-deg_CD45_GEX")
neoredp11 <- CreateSeuratObject(counts = neoredp11, min.features = 200, min.cells = 50)
neoredp11$patient <- "Patient6"
neoredp11$treatment <- "ADT+aCTLA4"
neoredp11$tissue <- "CD45+"
neoredp12 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/MC14_neo-RED-P_P005-C_CTLA4-deg_total_GEX")
neoredp12 <- CreateSeuratObject(counts = neoredp12, min.features = 200, min.cells = 50)
neoredp12$patient <- "Patient6"
neoredp12$treatment <- "ADT+aCTLA4"
neoredp12$tissue <- "Total"
neoredp13 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/MC17_neo-RED-P_P006-D_Deg_CD45_GEX")
neoredp13 <- CreateSeuratObject(counts = neoredp13, min.features = 200, min.cells = 50)
neoredp13$patient <- "Patient7"
neoredp13$treatment <- "ADT"
neoredp13$tissue <- "CD45+"
neoredp14 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/MC16_neo-RED-P_P006-D_Deg_total_GEX")
neoredp14 <- CreateSeuratObject(counts = neoredp14, min.features = 200, min.cells = 50)
neoredp14$patient <- "Patient7"
neoredp14$treatment <- "ADT"
neoredp14$tissue <- "Total"
neoredp15 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/MC19_neo-RED-P_P007_CD45_GEX")
neoredp15 <- CreateSeuratObject(counts = neoredp15, min.features = 200, min.cells = 50)
neoredp15$patient <- "Patient8"
neoredp15$treatment <- "ADT+aCTLA4"
neoredp15$tissue <- "CD45+"
neoredp16 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/MC18_neo-RED-P_P007_total_GEX")
neoredp16 <- CreateSeuratObject(counts = neoredp16, min.features = 200, min.cells = 50)
neoredp16$patient <- "Patient8"
neoredp16$treatment <- "ADT+aCTLA4"
neoredp16$tissue <- "Total"
neoredp17 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/MC22_neo-RED-P_P009-U_no-treatment_CD45_GEX")
neoredp17 <- CreateSeuratObject(counts = neoredp17, min.features = 200, min.cells = 50)
neoredp17$patient <- "Patient9"
neoredp17$treatment <- "Untreated"
neoredp17$tissue <- "CD45+"
neoredp18 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/MC21_neo-RED-P_P009-U_no-treatment_total_GEX")
neoredp18 <- CreateSeuratObject(counts = neoredp18, min.features = 200, min.cells = 50)
neoredp18$patient <- "Patient9"
neoredp18$treatment <- "Untreated"
neoredp18$tissue <- "Total"
neoredp19 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/MC24-neo-RED-P-P010-CD45-GEX")
neoredp19 <- CreateSeuratObject(counts = neoredp19, min.features = 200, min.cells = 50)
neoredp19$patient <- "Patient10"
neoredp19$treatment <- "ADT"
neoredp19$tissue <- "CD45+"
neoredp20 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/MC23-neo-RED-P-P010-total-GEX")
neoredp20 <- CreateSeuratObject(counts = neoredp20, min.features = 200, min.cells = 50)
neoredp20$patient <- "Patient10"
neoredp20$treatment <- "ADT"
neoredp20$tissue <- "Total"
# neoredp21<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/MC26-neo-RED-P-P011-CD45-GEX")
# neoredp21<-CreateSeuratObject(counts = neoredp21,min.features = 200,min.cells=50)
# neoredp21$patient="Patient11"
# neoredp21$treatment="Untreated"
# neoredp21$tissue="CD45+"
# neoredp22<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/MC25-neo-RED-P-P011-total-GEX")
# neoredp22<-CreateSeuratObject(counts = neoredp22,min.features = 200,min.cells=50)
# neoredp22$patient="Patient11"
# neoredp22$treatment="Untreated"
# neoredp22$tissue="Total"
neoredp23 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/MC20_neo-RED-P_P008_Untouched_GEX")
neoredp23 <- CreateSeuratObject(counts = neoredp23, min.features = 200, min.cells = 50)
neoredp23$patient <- "Patient12"
neoredp23$treatment <- "ADT+aCTLA4"
neoredp23$tissue <- "Total"
# neoredp24<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES/MC29_GEX/analysis/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES-MC29_GEX-cellranger-count-default/MC29_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
# neoredp24<-CreateSeuratObject(counts = neoredp24,min.features = 200,min.cells=50)
# neoredp24$patient="Patient13"
# neoredp24$treatment="Untreated"
# neoredp24$tissue="Total"
# neoredp25<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES/MC30_GEX/analysis/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES-MC30_GEX-cellranger-count-default/MC30_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
# neoredp25<-CreateSeuratObject(counts = neoredp25,min.features = 200,min.cells=50)
# neoredp25$patient="Patient13"
# neoredp25$treatment="Untreated"
# neoredp25$tissue="CD45+"
neoredp26 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES/MC31_GEX/analysis/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES-MC31_GEX-cellranger-count-default/MC31_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp26 <- CreateSeuratObject(counts = neoredp26, min.features = 200, min.cells = 50)
neoredp26$patient <- "Patient14"
neoredp26$treatment <- "ADT"
neoredp26$tissue <- "Total"
neoredp27 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES/MC32_GEX/analysis/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES-MC32_GEX-cellranger-count-default/MC32_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp27 <- CreateSeuratObject(counts = neoredp27, min.features = 200, min.cells = 50)
neoredp27$patient <- "Patient14"
neoredp27$treatment <- "ADT"
neoredp27$tissue <- "CD45+"
neoredp28 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES/MC33_GEX/analysis/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES-MC33_GEX-cellranger-count-default/MC33_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp28 <- CreateSeuratObject(counts = neoredp28, min.features = 200, min.cells = 50)
neoredp28$patient <- "Patient15"
neoredp28$treatment <- "ADT"
neoredp28$tissue <- "Total"
neoredp29 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES/MC34_GEX/analysis/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES-MC34_GEX-cellranger-count-default/MC34_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp29 <- CreateSeuratObject(counts = neoredp29, min.features = 200, min.cells = 50)
neoredp29$patient <- "Patient15"
neoredp29$treatment <- "ADT"
neoredp29$tissue <- "CD45+"
neoredp30 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES/MC35_GEX/analysis/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES-MC35_GEX-cellranger-count-default/MC35_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp30 <- CreateSeuratObject(counts = neoredp30, min.features = 200, min.cells = 50)
neoredp30$patient <- "Patient16"
neoredp30$treatment <- "ADT+aCTLA4"
neoredp30$tissue <- "Total"
neoredp31 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES/MC36_GEX/analysis/210913_BENJAMIN_WAN-I_2_HUMAN_10X_LANES-MC36_GEX-cellranger-count-default/MC36_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp31 <- CreateSeuratObject(counts = neoredp31, min.features = 200, min.cells = 50)
neoredp31$patient <- "Patient16"
neoredp31$treatment <- "ADT+aCTLA4"
neoredp31$tissue <- "CD45+"
neoredp32 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/210927_BENJAMIN_WAN-I_1_HUMAN_10X_LANE/MC37_GEX/analysis/210927_BENJAMIN_WAN-I_1_HUMAN_10X_LANE-MC37_GEX-cellranger-count-default/MC37_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp32 <- CreateSeuratObject(counts = neoredp32, min.features = 200, min.cells = 50)
neoredp32$patient <- "Patient17"
neoredp32$treatment <- "ADT"
neoredp32$tissue <- "Total"
neoredp33 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/210927_BENJAMIN_WAN-I_1_HUMAN_10X_LANE/MC38_GEX/analysis/210927_BENJAMIN_WAN-I_1_HUMAN_10X_LANE-MC38_GEX-cellranger-count-default/MC38_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp33 <- CreateSeuratObject(counts = neoredp33, min.features = 200, min.cells = 50)
neoredp33$patient <- "Patient17"
neoredp33$treatment <- "ADT"
neoredp33$tissue <- "CD45+"
neoredp34 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/220218_BENJAMIN_NEHA_2_HUMAN_10X/MC48-220124_P016-C-UNTOUCHED_GEX/MC48-220124_P016-C-UNTOUCHED_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp34 <- CreateSeuratObject(counts = neoredp34, min.features = 200, min.cells = 50)
neoredp34$patient <- "Patient18"
neoredp34$treatment <- "ADT+aCTLA4"
neoredp34$tissue <- "Total"
neoredp35 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/220218_BENJAMIN_NEHA_2_HUMAN_10X/MC49-220124-P016-C-CD45_GEX/MC49-220124-P016-C-CD45_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp35 <- CreateSeuratObject(counts = neoredp35, min.features = 200, min.cells = 50)
neoredp35$patient <- "Patient18"
neoredp35$treatment <- "ADT+aCTLA4"
neoredp35$tissue <- "CD45+"
neoredp36 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/220218_BENJAMIN_NEHA_2_HUMAN_10X/MC50-220131-P017-C-UNTOUCHED_GEX/MC50-220131-P017-C-UNTOUCHED_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp36 <- CreateSeuratObject(counts = neoredp36, min.features = 200, min.cells = 50)
neoredp36$patient <- "Patient19"
neoredp36$treatment <- "ADT+aCTLA4"
neoredp36$tissue <- "Total"
neoredp37 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/220218_BENJAMIN_NEHA_2_HUMAN_10X/MC51-220207-P01-C-UNTOUCHED_GEX/MC51-220207-P01-C-UNTOUCHED_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp37 <- CreateSeuratObject(counts = neoredp37, min.features = 200, min.cells = 50)
neoredp37$patient <- "Patient20"
neoredp37$treatment <- "ADT+aCTLA4"
neoredp37$tissue <- "Total"
neoredp38 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/220418_BENJAMIN_NEHA_3_HUMAN_10X/MC52-220224--P019-D-UNTOUCHED/analysis/220418_BENJAMIN_NEHA_3_HUMAN_10X-MC52-220224--P019-D-UNTOUCHED-cellranger-count-default/MC52-220224--P019-D-UNTOUCHED_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp38 <- CreateSeuratObject(counts = neoredp38, min.features = 200, min.cells = 50)
neoredp38$patient <- "Patient21"
neoredp38$treatment <- "ADT"
neoredp38$tissue <- "Total"
neoredp39 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/220418_BENJAMIN_NEHA_3_HUMAN_10X/MC53-220311--P020-D-UNTOUCHED/analysis/220418_BENJAMIN_NEHA_3_HUMAN_10X-MC53-220311--P020-D-UNTOUCHED-cellranger-count-default/MC53-220311--P020-D-UNTOUCHED_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp39 <- CreateSeuratObject(counts = neoredp39, min.features = 200, min.cells = 50)
neoredp39$patient <- "Patient22"
neoredp39$treatment <- "ADT"
neoredp39$tissue <- "Total"
# neoredp40=Read10X("/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/220513_BENJAMIN_NEHA_1_HUMAN_10X/MC58_220419_U006_U_CD45/MC58_220419_U006_U_CD45_cellranger_count_outs/filtered_feature_bc_matrix")
# neoredp40<-CreateSeuratObject(counts = neoredp40,min.features = 200,min.cells=50)
# neoredp40$patient="Patient23"
# neoredp40$treatment="Untreated"
# neoredp40$tissue="CD45+"
# neoredp41=Read10X("/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/220513_BENJAMIN_NEHA_1_HUMAN_10X/MC57_220419_U006_U_UNTOUCHED/MC57_220419_U006_U_UNTOUCHED_cellranger_count_outs/filtered_feature_bc_matrix")
# neoredp41<-CreateSeuratObject(counts = neoredp41,min.features = 200,min.cells=50)
# neoredp41$patient="Patient23"
# neoredp41$treatment="Untreated"
# neoredp41$tissue="Total"
# neoredp42=Read10X("/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/220901_BENJAMIN_NEHA_1_HUMAN_10X/MC59_GEX/analysis/220901_BENJAMIN_NEHA_1_HUMAN_10X-MC59_GEX-cellranger-count-default/MC59_GEX_cellranger_count_outs/filtered_feature_bc_matrix") #MC59_220802_U007_U_Untouched
# neoredp42<-CreateSeuratObject(counts = neoredp42,min.features = 200,min.cells=50)
# neoredp42$patient="Patient24"
# neoredp42$treatment="Untreated"
# neoredp42$tissue="Total"
# neoredp43=Read10X("/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/220901_BENJAMIN_NEHA_1_HUMAN_10X/MC60_GEX/analysis/220901_BENJAMIN_NEHA_1_HUMAN_10X-MC60_GEX-cellranger-count-default/MC60_GEX_cellranger_count_outs/filtered_feature_bc_matrix") #MC60_220802_U007_CD45
# neoredp43<-CreateSeuratObject(counts = neoredp43,min.features = 200,min.cells=50)
# neoredp43$patient="Patient24"
# neoredp43$treatment="Untreated"
# neoredp43$tissue="CD45+"
neoredp44 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/220901_BENJAMIN_NEHA_1_HUMAN_10X/MC61_GEX/analysis/220901_BENJAMIN_NEHA_1_HUMAN_10X-MC61_GEX-cellranger-count-default/MC61_GEX_cellranger_count_outs/filtered_feature_bc_matrix") # MC61_220808_P027_B_Untouched
neoredp44 <- CreateSeuratObject(counts = neoredp44, min.features = 200, min.cells = 50)
neoredp44$patient <- "Patient25"
neoredp44$treatment <- "ADT+aCTLA4"
neoredp44$tissue <- "Total"
# neoredp44<-Read10X(data.dir="/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/220901_BENJAMIN_NEHA_1_HUMAN_10X/MC62_GEX/analysis/220901_BENJAMIN_NEHA_1_HUMAN_10X-MC62_GEX-cellranger-count-default/MC62_GEX_cellranger_count_outs/filtered_feature_bc_matrix") #MC62_220812_U009_U_Untreated_Untouched
# neoredp44<-CreateSeuratObject(counts = neoredp44,min.features = 200,min.cells=50)
# neoredp44$patient="Patient26"
# neoredp44$treatment="Untreated"
# neoredp44$tissue="Total"
neoredp45 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/221005_BENJAMIN_NEHA_2_HUMAN_10X/MC63_GEX/analysis/221005_BENJAMIN_NEHA_2_HUMAN_10X-MC63_GEX-cellranger-count-default/MC63_GEX_cellranger_count_outs/filtered_feature_bc_matrix") # MC63_220902-PO26-A-ArmA_Untouched (degarelix only)
neoredp45 <- CreateSeuratObject(counts = neoredp45, min.features = 200, min.cells = 50)
neoredp45$patient <- "Patient27"
neoredp45$treatment <- "ADT"
neoredp45$tissue <- "Total"
neoredp46 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/221005_BENJAMIN_NEHA_2_HUMAN_10X/MC64_GEX/analysis/221005_BENJAMIN_NEHA_2_HUMAN_10X-MC64_GEX-cellranger-count-default/MC64_GEX_cellranger_count_outs/filtered_feature_bc_matrix") # MC64_NeoRed_control_090922022_U009-HA-Untouched
neoredp46 <- CreateSeuratObject(counts = neoredp46, min.features = 200, min.cells = 50)
neoredp46$patient <- "Patient28"
neoredp46$treatment <- "Untreated"
neoredp46$tissue <- "Total"
neoredp47 <- Read10X("/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/MC39/MC39_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp47 <- CreateSeuratObject(counts = neoredp47, min.features = 200, min.cells = 50)
neoredp47$patient <- "Patient29"
neoredp47$treatment <- "ADT+aCTLA4"
neoredp47$tissue <- "Total"
neoredp48 <- Read10X("/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/MC40/MC40_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp48 <- CreateSeuratObject(counts = neoredp48, min.features = 200, min.cells = 50)
neoredp48$patient <- "Patient29"
neoredp48$treatment <- "ADT+aCTLA4"
neoredp48$tissue <- "CD45+"
neoredp49 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/230127_BENJAMIN_PARIN_1_HUMAN_10X/MC79_GEX/analysis/230127_BENJAMIN_PARIN_1_HUMAN_10X-MC79_GEX-cellranger-count-default/MC79_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp49 <- CreateSeuratObject(counts = neoredp49, min.features = 200, min.cells = 50)
neoredp49$patient <- "Patient30"
neoredp49$treatment <- "Untreated"
neoredp49$tissue <- "Total"
neoredp50 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/230215_BENJAMIN_PARIN_1_HUMAN_SEQUENCING/MC80/analysis/230215_BENJAMIN_PARIN_1_HUMAN_SEQUENCING-MC80-cellranger-count-default/MC80_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp50 <- CreateSeuratObject(counts = neoredp50, min.features = 200, min.cells = 50)
neoredp50$patient <- "Patient31"
neoredp50$treatment <- "Untreated"
neoredp50$tissue <- "Total"
neoredp51 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/230215_BENJAMIN_PARIN_1_HUMAN_SEQUENCING/MC82/analysis/230215_BENJAMIN_PARIN_1_HUMAN_SEQUENCING-MC82-cellranger-count-default/MC82_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp51 <- CreateSeuratObject(counts = neoredp51, min.features = 200, min.cells = 50)
neoredp51$patient <- "Patient32"
neoredp51$treatment <- "Untreated"
neoredp51$tissue <- "Total"
neoredp52 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/230215_BENJAMIN_PARIN_1_HUMAN_SEQUENCING/MC83/analysis/230215_BENJAMIN_PARIN_1_HUMAN_SEQUENCING-MC83-cellranger-count-default/MC83_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp52 <- CreateSeuratObject(counts = neoredp52, min.features = 200, min.cells = 50)
neoredp52$patient <- "Patient33"
neoredp52$treatment <- "Untreated"
neoredp52$tissue <- "Total"
neoredp53 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/230215_BENJAMIN_PARIN_1_HUMAN_SEQUENCING/MC84_GEX/analysis/230215_BENJAMIN_PARIN_1_HUMAN_SEQUENCING-MC84_GEX-cellranger-count-default/MC84_GEX_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp53 <- CreateSeuratObject(counts = neoredp53, min.features = 200, min.cells = 50)
neoredp53$patient <- "Patient34"
neoredp53$treatment <- "Untreated"
neoredp53$tissue <- "Total"
neoredp54 <- Read10X(data.dir = "/Users/aleksandar/genomecenter.columbia.edu/ngs/release/singleCell/230215_BENJAMIN_PARIN_1_HUMAN_SEQUENCING/MC81/analysis/230215_BENJAMIN_PARIN_1_HUMAN_SEQUENCING-MC81-cellranger-count-default/MC81_cellranger_count_outs/filtered_feature_bc_matrix")
neoredp54 <- CreateSeuratObject(counts = neoredp54, min.features = 200, min.cells = 50)
neoredp54$patient <- "Patient35"
neoredp54$treatment <- "ADT+aCTLA4"
neoredp54$tissue <- "Total"


neoredp_list <- list(neoredp2, neoredp4, neoredp8, neoredp12, neoredp14, neoredp16, neoredp18, neoredp20, neoredp23, neoredp26, neoredp28, neoredp30, neoredp32, neoredp34, neoredp36, neoredp37, neoredp38, neoredp39, neoredp45, neoredp46, neoredp47, neoredp49, neoredp50, neoredp51, neoredp52, neoredp53, neoredp54)
rm(neoredp1, neoredp2, neoredp3, neoredp4, neoredp5, neoredp6, neoredp7, neoredp8, neoredp9, neoredp10, neoredp11, neoredp12, neoredp13, neoredp14, neoredp15, neoredp16, neoredp17, neoredp18, neoredp19, neoredp20, neoredp21, neoredp22, neoredp23, neoredp24, neoredp25, neoredp26, neoredp27, neoredp28, neoredp29, neoredp30, neoredp31, neoredp32, neoredp33, neoredp34, neoredp35, neoredp36, neoredp37, neoredp38, neoredp39, neoredp40, neoredp41, neoredp42, neoredp43, neoredp44, neoredp45, neoredp46, neoredp47, neoredp48, neoredp49, neoredp50, neoredp51, neoredp52, neoredp53, neoredp54)
for (i in 1:length(neoredp_list)) {
  p <- neoredp_list[[i]]
  p <- PercentageFeatureSet(p, pattern = "^MT-", col.name = "percent.mt")
  p <- subset(p, subset = percent.mt < 25 & nCount_RNA > 1500 & nCount_RNA < 10000) # >1000, <10000
  p <- SCTransform(p, vars.to.regress = c("nCount_RNA", "percent.mt"), return.only.var.genes = F, verbose = T, conserve.memory = T)
  p.singler <- CreateSinglerObject(p[["SCT"]]@counts,
    annot = NULL,
    project.name = "primecut", min.genes = 0,
    technology = "10X", species = "Human", citation = "",
    do.signatures = F, clusters = NULL, numCores = numCores,
    fine.tune = F, temp.dir = "/Users/aleksandar/Downloads",
    variable.genes = "de", reduce.file.size = T, do.main.types = T
  )
  p$hpca_labels <- p.singler$singler[[1]][[1]][[2]]
  p$hpca_main_labels <- p.singler$singler[[1]][[4]][[2]]
  p$blueprint_labels <- p.singler$singler[[2]][[1]][[2]]
  p$blueprint_main_labels <- p.singler$singler[[2]][[4]][[2]]
  p$hpca_pvals <- p.singler$singler[[1]][[1]][[3]]
  p$hpca_main_pvals <- p.singler$singler[[1]][[4]][[3]]
  p$blueprint_pvals <- p.singler$singler[[2]][[1]][[3]]
  p$blueprint_main_pvals <- p.singler$singler[[2]][[4]][[3]]
  neoredp_list[[i]] <- p
}
lapply(neoredp_list, ncol)
neoredp_list <- neoredp_list[which(lapply(neoredp_list, ncol) > 500)]
lapply(neoredp_list, ncol)
features <- SelectIntegrationFeatures(object.list = neoredp_list, nfeatures = 3000)
neoredp_list <- PrepSCTIntegration(object.list = neoredp_list, anchor.features = features, verbose = T)
neoredp_list <- lapply(X = neoredp_list, FUN = RunPCA, features = features)
anchors <- FindIntegrationAnchors(object.list = neoredp_list, normalization.method = "SCT", anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20, verbose = T, reference = 1)
rm(neoredp_list, features)
neoredp.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30, verbose = T)
rm(anchors)
neoredp.integrated$treatment <- factor(neoredp.integrated$treatment, levels = c("Untreated", "ADT", "ADT+aCTLA4"))

neoredp.integrated <- RunPCA(neoredp.integrated, features = VariableFeatures(object = neoredp.integrated))
neoredp.integrated <- RunUMAP(neoredp.integrated, dims = 1:50, verbose = FALSE, metric = "correlation")
neoredp.integrated <- FindNeighbors(neoredp.integrated, dims = 1:50, verbose = FALSE)
neoredp.integrated <- FindClusters(neoredp.integrated, resolution = seq(0.01, 1, by = 0.01), verbose = FALSE, algorithm = 1)
clust <- neoredp.integrated@meta.data[, which(grepl("integrated_snn_res.", colnames(neoredp.integrated@meta.data)))]
mat <- as.data.frame(t(neoredp.integrated$pca@cell.embeddings))
out <- sil_subsample(mat, clust)
means <- out[[1]]
sd <- out[[2]]
x <- seq(0.01, 1, by = 0.01)
errbar(x, means, means + sd, means - sd, ylab = "mean silhouette score", xlab = "resolution parameter")
lines(x, means)
best <- tail(x[which(means == max(means))], n = 1)
legend("topright", paste("Best", best, sep = " = "))
neoredp.integrated$seurat_clusters <- neoredp.integrated@meta.data[, which(colnames(neoredp.integrated@meta.data) == paste("integrated_snn_res.", best, sep = ""))]
Idents(neoredp.integrated) <- "seurat_clusters"
plot(DimPlot(neoredp.integrated, reduction = "umap", label = TRUE, label.size = 7, repel = T) + NoLegend())
markers <- FindAllMarkers(neoredp.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "wilcox")
top10 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
geneHeatmap_plot(neoredp.integrated@assays$SCT@scale.data, neoredp.integrated$seurat_clusters, top10$gene, n_top_genes_per_cluster = 5, scaled = F)
l <- neoredp.integrated$blueprint_labels
l[which(neoredp.integrated$blueprint_pvals > 0.1)] <- NA
l[which(l %in% names(which(table(l) < 150)))] <- NA
neoredp.integrated$l <- l
Idents(neoredp.integrated) <- "l"
plot(DimPlot(neoredp.integrated, reduction = "umap", label = TRUE, repel = T, label.size = 5) + NoLegend())
saveRDS(neoredp.integrated, file = "~/Documents/Documents/MED SCHOOL/summer 2019/neoredp.integrated.final.v3.rds")

## make metaCells -- target number of neighbors should be 10000/(mean UMI count) rounded to the nearest integer. Will vary by depth/quality of the dataset. If number of neighbors is too high be cautious of cell-type-mixing and exclude small cell clusters by adjusting sizeThresh.
neoredp.integrated.meta <- MakeCMfA(dat.mat = as.matrix(neoredp.integrated[["SCT"]]@counts), clustering = neoredp.integrated$seurat_clusters, out.dir = "/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/prime-cut-metacells/", out.name = "neoredp_integrated", sizeThresh = 50, numNeighbors = 5)

### THIS IS WHERE YOU RUN ARACNe ON THE METACELL OUTPUT

## run VIPER
## load aracne networks
filenames <- list.files("/Users/aleksandar/Documents/Documents/MSK prostate data/single_cell_prostate_nets", pattern = "*.rds", full.names = TRUE)
nets <- lapply(filenames, readRDS)
dat <- neoredp.integrated@assays$integrated@scale.data
### chunked in groups of 400 cells at a time to avoid overloading memory
indices <- seq(1, ncol(dat), by = 400)
seurat_viper_list <- list()
for (i in 1:(length(indices) - 1)) {
  vp <- viper(dat[, indices[i]:indices[i + 1]], nets, method = "none")
  seurat_viper_list <- c(seurat_viper_list, list(vp))
}
vp <- viper(dat[, indices[i + 1]:ncol(dat)], nets, method = "none")
seurat_viper_list <- c(seurat_viper_list, list(vp))
rm(dat)
vp_meta <- seurat_viper_list[[1]]
for (i in 2:length(seurat_viper_list)) {
  vp_meta <- cbind(vp_meta, seurat_viper_list[[i]])
}
vp_meta <- vp_meta[, colnames(neoredp.integrated)]


## VIPER clustering
neoredp.integrated.meta.vp <- CreateSeuratObject(counts = vp_meta[cbcMRs, ])
neoredp.integrated.meta.vp@assays$RNA@scale.data <- as.matrix(neoredp.integrated.meta.vp@assays$RNA@data)
neoredp.integrated.meta.vp <- RunPCA(neoredp.integrated.meta.vp, features = rownames(neoredp.integrated.meta.vp))
neoredp.integrated.meta.vp <- RunUMAP(neoredp.integrated.meta.vp, dims = 1:30, verbose = FALSE, metric = "correlation")
neoredp.integrated.meta.vp <- FindNeighbors(neoredp.integrated.meta.vp, dims = 1:30, verbose = FALSE)
neoredp.integrated.meta.vp <- FindClusters(neoredp.integrated.meta.vp, resolution = seq(0.01, 1, by = 0.01), verbose = FALSE, algorithm = 1)
clust <- neoredp.integrated.meta.vp@meta.data[, which(grepl("RNA_snn_res.", colnames(neoredp.integrated.meta.vp@meta.data)))]
mat <- as.data.frame(t(neoredp.integrated.meta.vp$pca@cell.embeddings))
out <- sil_subsample(mat, clust)
means <- out[[1]]
sd <- out[[2]]
x <- seq(0.01, 1, by = 0.01)
errbar(x, means, means + sd, means - sd, ylab = "mean silhouette score", xlab = "resolution parameter")
lines(x, means)
best <- tail(x[which(means[15:length(means)] == max(means[15:length(means)])) + 14], n = 1)
legend("topright", paste("Best", best, sep = " = "))
neoredp.integrated.meta.vp$seurat_clusters <- neoredp.integrated.meta.vp@meta.data[, which(colnames(neoredp.integrated.meta.vp@meta.data) == paste("RNA_snn_res.", best, sep = ""))]
Idents(neoredp.integrated.meta.vp) <- "seurat_clusters"
neoredp.integrated.meta.vp$patient <- neoredp.integrated$patient[colnames(neoredp.integrated.meta.vp)]
neoredp.integrated.meta.vp$treatment <- neoredp.integrated$treatment[colnames(neoredp.integrated.meta.vp)]
neoredp.integrated.meta.vp$tissue <- neoredp.integrated$tissue[colnames(neoredp.integrated.meta.vp)]
neoredp.integrated.meta.vp$blueprint_labels <- neoredp.integrated$blueprint_labels[colnames(neoredp.integrated.meta.vp)]
neoredp.integrated.meta.vp$blueprint_pvals <- neoredp.integrated$blueprint_pvals[colnames(neoredp.integrated.meta.vp)]
neoredp.integrated.meta.vp$l <- neoredp.integrated$l[colnames(neoredp.integrated.meta.vp)]
neoredp.integrated.meta.vp$treatment <- factor(neoredp.integrated.meta.vp$treatment, levels = c("Untreated", "ADT", "ADT+aCTLA4"))
# neoredp.integrated.meta.vp$seurat_clusters=mapvalues(neoredp.integrated.meta.vp$seurat_clusters, from = 1:length(unique(neoredp.integrated.meta.vp$seurat_clusters))-1, to = c("T-cell","Tumor.1","Tumor.2","Endothelial","Myeloid","B-cell","Adipocyte","Fibroblast.1","Fibroblast.2","NK-cell","Tumor.3","MEP"))
Idents(neoredp.integrated.meta.vp) <- "seurat_clusters"
plot(DimPlot(neoredp.integrated.meta.vp, reduction = "umap", label = TRUE, label.size = 7, repel = T) + NoLegend())
# model <- lda(t(as.matrix(neoredp.integrated.meta.vp@assays$RNA@counts)),neoredp.integrated.meta.vp$seurat_clusters)
# lda=predict(model)$x
# neoredp.integrated.meta.vp@reductions[["lda"]]=CreateDimReducObject(embeddings = lda, key = "LD_", assay = DefaultAssay(neoredp.integrated.meta.vp))
# DimPlot(neoredp.integrated.meta.vp,reduction = "lda",label = T,repel=T)
markers.vp <- FindAllMarkers(neoredp.integrated.meta.vp, only.pos = TRUE, min.pct = 0, logfc.threshold = 0, test.use = "t")
top10 <- markers.vp %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
geneHeatmap_plot(neoredp.integrated.meta.vp@assays$RNA@counts, neoredp.integrated.meta.vp$seurat_clusters, top10$gene, n_top_genes_per_cluster = 5, scaled = T)
l <- neoredp.integrated.meta.vp$blueprint_labels
l[which(neoredp.integrated.meta.vp$blueprint_pvals > 0.1)] <- NA
l[which(l %in% names(which(table(l) < 150)))] <- NA
neoredp.integrated.meta.vp$l <- l
plot(DimPlot(neoredp.integrated.meta.vp, reduction = "umap", label = TRUE, repel = T, label.size = 5, group.by = "l") + NoLegend())
saveRDS(neoredp.integrated.meta.vp, file = "~/Documents/Documents/MED SCHOOL/summer 2019/neoredp.integrated.meta.vp.final.v3.rds")
