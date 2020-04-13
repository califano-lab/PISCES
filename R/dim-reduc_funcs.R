#' Builds a UMAP transformation based on the shared, cell-by-cell master regulators. Utilizes custom UMAP parameters ('pearson', 25 neighbors).
#'
#' @param dat.mat Matrix of protein activity (features X samples).
#' @param num.mrs Number of top master regulators to take from each cell. Default of 50.
#' @return Returns a UMAP based on the cell-by-cell master regulators.
#' @export
cbcMR_UMAP <- function(dat.mat, num.mrs = 50) {
  # set UMAP parameters
  umap_custom <- umap::umap.defaults
  umap_custom$n_neighbors <- 25
  umap_custom$metric <- 'pearson'
  # identify MRs
  cbc.mrs <- apply(dat.mat, 2, function(x) { names(sort(x, decreasing = TRUE))[1:num.mrs] })
  cbc.mrs <- unique(unlist(as.list(cbc.mrs)))
  # generate and return UMAP
  cbc.umap <- umap::umap(t(dat.mat[ match(cbc.mrs, rownames(dat.mat)) , ]), config = umap_custom, init = "random")
  return(cbc.umap)
}

#' Generates a UMAP based on a set of proteins.
#'
#' @param dat.mat Matrix of protein activity (features X samples).
#' @param mrs List of proteins to use as master regulators.
#' @return UMAP object.
#' @export
CustomUMAP <- function(dat.mat) {
  require(umap)
  # set UMAP parameters
  set.seed(1)
  umap_custom <- umap::umap.defaults
  umap_custom$n_neighbors <- 25
  umap_custom$metric <- 'pearson'
  # compute umap
  c.umap <- umap::umap(t(dat.mat), config = umap_custom, init = 'random')
  return(c.umap)
}
