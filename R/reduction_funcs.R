#' Generates the specified number of PCs using RSpectra.
#' 
#' @param dat.mat Matrix of data (features X samples).
#' @param num.pcs Number of PCs to be generated. Default of 10.
#' @return A list, with with pca matrix (x), a vector of feature sds (sdev), and the rotation matrix (rotation).
#' @export
fast_pca <- function(dat.mat, num.pcs = 10) {
  require(RSpectra)
  # generate svd
  dat.mat <- t(dat.mat)
  svd.obj <- svds(dat.mat, k = min(num.pcs, ncol(dat.mat)), nu = 0,
                  opts = list(center = TRUE, scale = TRUE))
  # format pca matrix
  pca.mat <- dat.mat %*% svd.obj$v
  rownames(pca.mat) <- rownames(dat.mat)
  colnames(pca.mat) <- paste('PC', 1:num.pcs, sep = '.')
  # return pca list
  pca.list <- list('x' = pca.mat, 'sdev' = svd.obj$d, 'rotation' = svd.obj$v)
  return(pca.list)
}

#' Generates a UMAP using the uwot package.
#' 
#' @param dat.mat Matrix of data (features X samples).
#' @param dense.map Flag to use dense UMAP. Default of TRUE.
#' @return Two column matrix with UMAP representation of samples.
#' @export
make_umap <- function(dat.mat, dense = TRUE) {
  require(uwot)
  require(densvis)
  if (dense) {
    umap.mat <- densmap(t(dat.mat))
  } else {
    umap.mat <- umap(t(dat.mat), metric = 'correlation')
  }
  rownames(umap.mat) <- colnames(dat.mat)
  return(umap.mat)
}

#' Generates an MDS using base R functionality.
#' 
#' @param dist.mat Distance matrix.
#' @return Two column MDS matrix.
#' @export
make_mds <- function(dist.mat) {
  mds.mat <- cmdscale(dist.mat, k = 2)
  return(mds.mat)
}

#' Generates a distance matrix using the square root of (1 - spearman).
#' 
#' @param dat.mat Matrix of data (features X samples).
#' @return Distance structure.
#' @export
cor_dist <- function(dat.mat) {
  dist.mat <- as.dist(sqrt(1 - cor(dat.mat, method = 'spearman')))
  return(dist.mat)
}
