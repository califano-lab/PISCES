#' Generates a KNN matrix.
#'
#' @param dist.mat Distance matrix.
#' @param k Number of neighbors. Default of 5.
#' @return A naighbor matrix with samples in rows and neighbors in columns.
#' @export
KNN <- function(dist.mat, k = 5){
  dist.mat <- as.matrix(dist.mat)
  n <- nrow(dist.mat)
  neighbor.mat <- matrix(0L, nrow = n, ncol = k)
  for (i in 1:n) {
    neighbor.mat[i,] <- order(dist.mat[i,])[2:(k + 1)]
  }
  return(neighbor.mat)
}

#' Generates meta cells for each cluster
#' 
#' @param counts.mat Matrix of filtered counts (genes X samples).
#' @param dist.mat Distance matrix for this data.
#' @param clust.vect Clustering vector. If not specified, entire matrix will be used.
#' @param num.neighbors Number of neighbors to use in metacells. Default of 5.
#' @param subset Switch to control subsetting to 250 samples for ARACNe. Default of TRUE.
#' @param cpm.norm Switch to control cpm normalization of meta cell matrix. Default of TRUE.
#' @param min.samps Minimum number of samples in a cluster required for meta cells. Default of 500. 
#' @return A list of meta cell matrices for all clusters with enough samples.
#' @export
MetaCells <- function(counts.mat, dist.mat, clust.vect, num.neighbors = 5, subset = TRUE, cpm.norm = TRUE, min.samps = 500) {
  dist.mat <- as.matrix(dist.mat)
  # make dummy clusering if missing
  if (missing(clust.vect)) {
    clust.vect <- rep(1, ncol(counts.mat))
  }
  # make meta cell matrix for each cluster
  clust.table <- table(clust.vect)
  meta.mats <- list()
  for (clust.name in names(clust.table)) {
    if (clust.table[clust.name] > min.samps) {
      # get cluster matrix
      clust.samps <- which(clust.vect == clust.name)
      clust.mat <- counts.mat[, clust.samps]
      clust.dist <- as.dist(dist.mat[clust.samps, clust.samps])
      # get nearest neighbors
      knn.neighbors <- KNN(clust.dist, k = num.neighbors)
      # impute
      imp.mat <- matrix(0, nrow = nrow(clust.mat), ncol = ncol(clust.mat))
      rownames(imp.mat) <- rownames(clust.mat); colnames(imp.mat) <- colnames(clust.mat)
      for (j in 1:ncol(clust.mat)) {
        neighbor.mat <- clust.mat[, c(j, knn.neighbors[j,])]
        imp.mat[,j] <- rowSums(neighbor.mat)
      }
      # subset, cpm, and store
      if (subset) { imp.mat <- imp.mat[, sample(colnames(imp.mat), 250)] }
      if (cpm.norm) { imp.mat <- CPMTransform(imp.mat) }
      meta.mats[[clust.name]] <- imp.mat
    }
  }
  return(meta.mats)
}