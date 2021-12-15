#' Adds a named assay 'PISCES' to a Seurat object.
#' 
#' @param seurat.obj Seurat object w/ RNA assay.
#' @return Seurat object w/ original assay and added 'PISCES' assay.
#' @export
AddPISCESAssay <- function(seurat.obj) {
  pisces.assay <- CreateAssayObject(counts = seurat.obj@assays$RNA@counts)
  seurat.obj[['PISCES']] <- pisces.assay
  seurat.obj@active.assay <- 'PISCES'
  return(seurat.obj)
}

#' Checks if the PISCES assay is present in a given object.
#' 
#' @param data.object Seurat object.
#' @return Boolean TRUE or FALSE.
#' @export
HasPISCESAssay <- function(data.object) {
  return('PISCES' %in% names(data.object@assays))
}

#' Performas a CPM normalization on the counts in the given seurat object or matrix.
#' If a seurat object, stores results in the data field; otherwise, returns cpm matrix.
#'
#' @param data.object Either a Seurat object or a matrix of raw count data (genes X samples).
#' @param l2 Optional log2 normalization switch. Default of False.
#' @param remove.zeros Removes genes with zero expression across all samples. True by default.
#' @param pseudo Optional pseudo count logical. Default of False.
#' @return CPM matrix, or appropriately adjusted seurat object.
#' @export
CPMTransform <- function(data.object, remove.zeros = TRUE, l2 = FALSE, pseudo = FALSE) {
  # check if seurat object
  if (class(data.object)[1] == "Seurat") {
    if (!HasPISCESAssay(data.object)) {
      AddPISCESAssay(data.object)
    }
    dat.mat <- as.matrix(data.object@assays$PISCES@counts)
  } else {
    dat.mat <- data.object
  }
  # pseudo count if specified
  if (pseudo) { dat.mat <- dat.mat + 1 }
  # remove genes with zero reads
  if (remove.zeros) { dat.mat <- dat.mat[which(rowSums(dat.mat) > 0),] }
  # cpm transform
  cpm.mat <- t(t(dat.mat) / (colSums(dat.mat) / 1e6))
  # log2 if specified
  if (l2) { cpm.mat <- log2(cpm.mat + 1) } 
  # return or add to object
  if (class(data.object)[1] == "Seurat") {
    data.object@assays$PISCES@data <- cpm.mat 
    return(data.object)
  } else {
    return(cpm.mat)
  }
}

#' Generates a gene expression signature (GES) using internal normalization.
#' If a seurat object, stores results in misc$GES; otherwise, returns GES matrix.
#'
#' @param data.object Either a Seurat object or a matrix of raw count data (genes X samples).
#' @return GES matrix, or appropriately adjusted seurat object.
#' @export
GESTransform <- function(data.object) {
  # check if seurat object
  if (class(data.object)[1] == "Seurat") {
    if (!HasPISCESAssay(data.object)) {
      AddPISCESAssay(data.object)
    }
    dat.mat <- as.matrix(data.object@assays$PISCES@data)
  } else {
    dat.mat <- data.object
  }
  # generate GES
  ges.mat <- t(apply(dat.mat, 1, function(x) {
    (x - mean(x)) / sd(x)
  }))
  # return
  if (class(data.object)[1] == "Seurat") {
    data.object@assays$PISCES@misc[['GES']] <- ges.mat 
    return(data.object)
  } else {
    return(ges.mat)
  }
}

#' Performs a rank transformation on a given matrix, typically as an alternative GES generation technique.
#' If a seurat object, stores results in misc$GES; otherwise, returns GES matrix.
#'
#' @param data.object Either a Seurat object or a matrix of raw count data (genes X samples).
#' @return Rank transformed matrix, or appropriately adjusted seurat object.
#' @export
RankGES <- function(data.object) {
  # check if seurat object
  if (class(data.object)[1] == "Seurat") {
    if (!HasPISCESAssay(data.object)) {
      AddPISCESAssay(data.object)
    }
    dat.mat <- as.matrix(data.object@assays$PISCES@counts)
  } else {
    dat.mat <- data.object
  }
  # generate transformation
  rank.mat <- apply(dat.mat, 2, rank)
  median <- apply(rank.mat, 1, median)
  mad <- apply(rank.mat, 1, mad)
  rank.mat <- (rank.mat - median) / mad
  # return
  if (class(data.object)[1] == "Seurat") {
    data.object@assays$PISCES@misc[['GES']] <- rank.mat 
    return(data.object)
  } else {
    return(rank.mat)
  }
}

#' Generates a distance matrix using sqrt(1-cor(x)) as a distance metric.
#' If a seurat object, stores results in cor.dist (misc); otherwise, returns distance matrix.
#' 
#' @param data.object Either a Seurat object or a matrix of data (features X samples).
#' @param cor.method Method argument for cor function; spearman by default.
#' @return Rank transformed matrix, or appropriately adjusted seurat object.
#' @export
CorDist <- function(data.object, cor.method = 'spearman') {
  # check if seurat object
  if (class(data.object)[1] == "Seurat") {
    dat.mat <- as.matrix(data.object@assays[[data.object@active.assay]]@scale.data)
  } else {
    dat.mat <- data.object
  }
  # generate distane matrix
  dist.mat <- as.dist(sqrt(1 - cor(dat.mat, method = cor.method)))
  # return
  if (class(data.object)[1] == "Seurat") {
    data.object@assays[[data.object@active.assay]]@misc[['dist.mat']] <- dist.mat 
    return(data.object)
  } else {
    return(dist.mat)
  }
}