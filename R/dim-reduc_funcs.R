#' Builds a UMAP using the uwot package. Optionally can subset to most significant data within each sample.
#' 
#' @param dat.mat Matrix of data (features X samples).
#' @param cbc Flag to perform cell-by-cell UMAP. Default of FALSE.
#' @param num.feats Number of features to use if cbc flag is TRUE. Default of 10.
#' @return A 2-column mantrix w/ UMAP coordinates.
#' @export
MakeUMAP <- function(dat.mat, cbc = FALSE, num.feats = 10) {
  # adjust data if cbc is specified
  if (cbc) {
    cbc.feats <- apply(dat.mat, 2, function(x) { names(sort(x, decreasing = TRUE))[1:num.feats] })
    cbc.feats <- unique(unlist(as.list(cbc.feats)))
    dat.mat <- dat.mat[match(cbc.feats, rownames(dat.mat)) ,]
  }
  # make umap
  umap.mat <- uwot::umap(t(dat.mat), metric = 'correlation')
  rownames(umap.mat) <- colnames(dat.mat)
  return(umap.mat)
}
