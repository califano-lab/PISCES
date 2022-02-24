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
  rownames(neighbor.mat) <- colnames(dist.mat)
  return(neighbor.mat)
}

#' Generates meta cells for each cluster
#' 
#' @param counts.mat Matrix of filtered counts (genes X samples).
#' @param dist.mat Distance matrix for this data.
#' @param clust.vect Clustering vector. If not specified, entire matrix will be used.
#' @param num.neighbors Number of neighbors to use in metacells. Default of 5.
#' @param subset Number of cells to subset to. Default of 250. No subsetting if set equal to NULL.
#' @param min.samps Minimum number of samples in a cluster required for meta cells. Default of 500. 
#' @return A list of meta cell matrices for all clusters with enough samples.
#' @export
MetaCells <- function(counts.mat, dist.mat, clust.vect, num.neighbors = 5, subset = 250, min.samps = 500) {
  n <- subset
  counts.mat <- as.matrix(counts.mat)
  dist.mat <- as.matrix(dist.mat)
  # dummy clustering vector if not specified
  if (missing(clust.vect)) {
    clust.vect <- rep(1, ncol(counts.mat))
    names(clust.vect) <- colnames(counts.mat)
  }
  clust.labels <- sort(unique(clust.vect))
  clust.labels <- as.character(clust.labels)
  # metacell matrix for each cluster
  meta.mats <- list()
  for (cl in clust.labels) {
    clust.samps <- names(clust.vect)[which(clust.vect == cl)]
    if (length(clust.samps) > min.samps) {
      print(paste("Making metacell matrix for cluster ", cl, "...", sep = ''))
      # get cluster objects
      clust.counts <- counts.mat[,clust.samps]
      clust.dist <- dist.mat[clust.samps, clust.samps]
      knn.mat <- KNN(clust.dist, k = num.neighbors)
      if (is.null(subset)) {
        sub.samps <- clust.samps
        n <- length(sub.samps)
      } else {
        sub.samps <- sample(clust.samps, subset)
      }
      # impute matrix
      imp.mat <- matrix(0L, nrow = nrow(clust.counts), ncol = n)
      rownames(imp.mat) <- rownames(counts.mat); colnames(imp.mat) <- sub.samps
      for (ss in sub.samps) {
        neighbor.vect <- c(ss, rownames(knn.mat)[knn.mat[ss,]])
        ss.mat <- clust.counts[, neighbor.vect]
        imp.mat[,ss] <- rowSums(ss.mat)
      }
      # normalize
      imp.mat <- CPMTransform(imp.mat)
      meta.mats[[as.character(cl)]] <- imp.mat
    }
  }
  return(meta.mats)
}
