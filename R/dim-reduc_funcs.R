#' Builds a UMAP using the uwot package. Optionally can subset to most significant data within each sample.
#' 
#' @param dat.object Seurat object w/ PISCES assay and viper matrix in 'misc'.
#' @param cbc Flag to perform cell-by-cell UMAP. Default of FALSE.
#' @param num.feats Number of features to use if cbc flag is TRUE. Default of 10.
#' @return Modified Seurat object w/ UMAP object in misc$umap
#' @export
MakeUMAP <- function(dat.object, cbc = FALSE, num.feats = 10) {
  # check if seurat object
  if (class(dat.object)[1] == "Seurat") {
    dat.mat <- dat.object@assays[[dat.object@active.assay]]@scale.data
  } else {
    dat.mat <- dat.object
  }
  # adjust data if cbc is specified
  if (cbc) {
    cbc.feats <- apply(dat.mat, 2, function(x) { names(sort(x, decreasing = TRUE))[1:num.feats] })
    cbc.feats <- unique(unlist(as.list(cbc.feats)))
    dat.mat <- dat.mat[match(cbc.feats, rownames(dat.mat)) ,]
  }
  # make umap
  umap.mat <- uwot::umap(t(dat.mat), metric = 'correlation')
  rownames(umap.mat) <- colnames(dat.mat)
  # check if seurat object
  if (class(dat.object)[1] == "Seurat") {
    dat.object@assays[[dat.object@active.assay]]@misc[['umap']] <- umap.mat
    return(dat.object)
  } else {
    return(umap.mat)
  }
}
#' Builds a UMAP using the uwot package. Optionally can subset to most significant data within each sample.
#' 
#' @param dat.object Seurat object w/ PISCES assay and viper matrix in 'misc'.
#' @param num.dims Number of dimensions to return. Default of 2.
#' @return Modified Seurat object w/ UMAP object in misc$umap
#' @export
MakeMDS <- function(dat.object, num.dims = 2) {
  # check if seurat object
  if (class(dat.object)[1] == "Seurat") {
    dist.mat <- dat.object@assays[[dat.object@active.assay]]@misc[['dist.mat']]
  } else {
    dist.mat <- dat.object
  }
  # generate mds
  mds.mat <- cmdscale(dist.mat, k = num.dims)
  # check if seurat object
  if (class(dat.object)[1] == "Seurat") {
    dat.object@assays[[dat.object@active.assay]]@misc[['mds']] <- mds.mat
    return(dat.object)
  } else {
    return(mds.mat)
  }
}

#' Gets the number of PCA features to use from a Seurat object based on the given variance theshold.
#' 
#' @param seurat.obj Seurat object with PCA reduction.
#' @param var.thresh Amount of variance to include. Default of 0.9.
#' @return Number of PCA features to use.
#' @export
GetPCAFeats <- function(seurat.obj, var.thresh = 0.9) {
  pca.var <- seurat.obj@reductions$pca@stdev**2; pca.var <- pca.var / sum(pca.var)
  var.sum <- cumsum(pca.var)
  num.dims <- tail(which(var.sum < var.thresh), 1) + 1
  return(num.dims)
}