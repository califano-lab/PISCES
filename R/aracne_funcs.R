#' Generates meta cells for each cluster
#' 
#' @param counts.mat Matrix of filtered counts (genes X samples).
#' @param dist.mat Distance matrix for this data.
#' @param clust.vec Clustering vector. If not specified, entire matrix will be used.
#' @param num.neighbors Number of neighbors to use in metacells. Default of 5.
#' @param subset Number of cells to subset to. Default of 250. No subsetting if set equal to NULL.
#' @param min.samps Minimum number of samples in a cluster required for meta cells. Default of 500. 
#' @return A list of meta cell matrices for all clusters with enough samples.
#' @export
make_metacells <- function(counts.mat, dist.mat, clust.vec, num.neighbors = 5, subset = 250, min.samps = 500) {
  # ensure matrix format
  counts.mat <- as.matrix(counts.mat)
  dist.mat <- as.matrix(dist.mat)
  # dummy clustering vector if not specified
  if (missing(clust.vec)) {
    clust.vec <- rep(1, ncol(counts.mat))
    names(clust.vec) <- colnames(counts.mat)
  }
  clust.labels <- sort(unique(clust.vec))
  clust.labels <- as.character(clust.labels)
  # make metacell matrix for each cluster
  meta.mats <- list()
  for (cl in clust.labels) {
    clust.samps <- names(clust.vec)[which(clust.vec == cl)]
    if (length(clust.samps) > min.samps) {
      print(paste("Making metacell matrix for cluster ", cl, "...", sep = ''))
      # get cluster objects
      clust.counts <- counts.mat[,clust.samps]
      clust.dist <- dist.mat[clust.samps, clust.samps]
      knn.mat <- knn(clust.dist, k = num.neighbors)
      if (is.null(subset)) {
        sub.samps <- clust.samps
        sub.size <- length(sub.samps)
      } else {
        sub.size <- min(length(clust.samps), subset)
        sub.samps <- sample(clust.samps, sub.size)
      }
      # impute matrix
      imp.mat <- matrix(0L, nrow = nrow(clust.counts), ncol = sub.size)
      rownames(imp.mat) <- rownames(counts.mat); colnames(imp.mat) <- sub.samps
      for (ss in sub.samps) {
        neighbor.vect <- c(ss, rownames(knn.mat)[knn.mat[ss,]])
        ss.mat <- clust.counts[, neighbor.vect]
        imp.mat[,ss] <- rowSums(ss.mat)
      }
      print(dim(imp.mat))
      # normalize
      imp.mat <- pflpf_norm(imp.mat)
      meta.mats[[as.character(cl)]] <- imp.mat
    }
  }
  return(meta.mats)
}

#' Generates a KNN matrix.
#'
#' @param dist.mat Distance matrix.
#' @param k Number of neighbors. Default of 5.
#' @return A naighbor matrix with samples in rows and neighbors in columns.
#' @export
knn <- function(dist.mat, k = 5){
  dist.mat <- as.matrix(dist.mat)
  n <- nrow(dist.mat)
  neighbor.mat <- matrix(0L, nrow = n, ncol = k)
  for (i in 1:n) {
    neighbor.mat[i,] <- order(dist.mat[i,])[2:(k + 1)]
  }
  rownames(neighbor.mat) <- colnames(dist.mat)
  return(neighbor.mat)
}

#' Saves a matrix in a format for input to ARACNe
#'
#' @param dat.mat Matrix of data (genes X samples).
#' @param out.file Output file where matrix will be saved.
#' @param subset Switch for subsetting the matrix to 500 samples. Default TRUE.
#' @export
ARACNeTable <- function(dat.mat, out.file, subset = TRUE) {
  dat.mat <- dat.mat[!duplicated(rownames(dat.mat)), ]
  if (subset) {
    dat.mat <- dat.mat[, sample(colnames(dat.mat), min(ncol(dat.mat), 500)) ]
  }
  sample.names <- colnames(dat.mat)
  gene.ids <- rownames(dat.mat)
  m <- dat.mat
  mm <- rbind( c("gene", sample.names), cbind(gene.ids, m))
  write.table( x = mm , file = out.file ,
               sep="\t", quote = F , row.names = F , col.names = F )
}