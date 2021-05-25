#' Generates basic quality control plots from raw gene expression data.
#'
#' @param seurat.obj A Seurat object w/ MT% added to metadata.
#' @param plot.path Optional argumetn of save path for plot.
#' @export
QCPlots <- function(seurat.obj, plot.path) {
  samp.factor <- as.factor(rep('raw', ncol(seurat.obj)))
  ## sequencing depth plot
  p1.dat <- data.frame('Depth' = seurat.obj$nCount_RNA, 'Sample' = samp.factor)
  p1 <- ggplot2::ggplot(p1.dat, ggplot2::aes(x=Sample, y=Depth)) + ggplot2::geom_violin(color = '#F8766D', fill = '#F8766D') +
    ggplot2::ylab('Depth') + ggplot2::xlab('') + ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x=ggplot2::element_blank())
  ## detected gene plot
  p2.dat <- data.frame('dgenes' = seurat.obj$nFeature_RNA, 'Sample' = samp.factor)
  p2 <- ggplot2::ggplot(p2.dat, ggplot2::aes(x=Sample, y=dgenes)) + ggplot2::geom_violin(color = '#00BA38', fill = '#00BA38') +
    ggplot2::ylab('Detected Genes') + ggplot2::xlab('') + ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x=ggplot2::element_blank())
  ## mt percentage plot
  p3.dat <- data.frame('mt' = seurat.obj$percent.mt, 'Sample' = samp.factor)
  p3 <- ggplot2::ggplot(p3.dat, ggplot2::aes(x=Sample, y=mt)) + ggplot2::geom_violin(color = '#619CFF', fill = '#619CFF') +
    ggplot2::ylab('MT%') + ggplot2::xlab('') + ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x=ggplot2::element_blank())
  ## arrange and plot
  if (!missing(plot.path)) {
    jpeg(plot.path, height = 600, width = 1000)
    print(ggpubr::ggarrange(plotlist = list(p1, p2, p3), ncol = 3))
    dev.off()
  } else {
    ggpubr::ggarrange(plotlist = list(p1, p2, p3), ncol = 3)
  }
}

#' Performas a CPM normalization on the counts in the given seurat object or matrix.
#' If a seurat object, stores results in the data field; otherwise, returns cpm matrix.
#'
#' @param data.object Either a Seurat object or a matrix of raw count data (genes X samples).
#' @param l2 Optional log2 normalization switch. Default of False.
#' @param pseudo Optional pseudo count logical. Default of False.
#' @return CPM matrix, or appropriately adjusted seurat object.
#' @export
CPMTransform <- function(data.object, l2 = FALSE, pseudo = FALSE) {
  # check if seurat object
  if (class(data.object)[1] == "Seurat") {
    dat.mat <- as.matrix(seurat.obj@assays$RNA@counts)
  } else {
    dat.mat <- data.object
  }
  # pseudo count if specified
  if (pseudo) { dat.mat <- dat.mat + 1 }
  # cpm transform
  cpm.mat <- t(t(dat.mat) / (colSums(dat.mat) / 1e6))
  # log2 if specified
  if (l2) { cpm.mat <- log2(cpm.mat + 1) } 
  # return or add to object
  if (class(data.object)[1] == "Seurat") {
    seurat.obj@assays$RNA@data <- cpm.mat 
    return(seurat.obj)
  } else {
    return(cpm.mat)
  }
}

#' Generates a gene expression signature (GES) using internal normalization.
#' If a seurat object, stores results in scale.data field; otherwise, returns GES matrix.
#'
#' @param data.object Either a Seurat object or a matrix of raw count data (genes X samples).
#' @return GES matrix, or appropriately adjusted seurat object.
#' @export
GESTransform <- function(data.object) {
  # check if seurat object
  if (class(data.object)[1] == "Seurat") {
    dat.mat <- as.matrix(seurat.obj@assays$RNA@data)
  } else {
    dat.mat <- data.object
  }
  # generate GES
  ges.mat <- t(apply(dat.mat, 1, function(x) {
    (x - mean(x)) / sd(x)
  }))
  # return
  if (class(data.object)[1] == "Seurat") {
    seurat.obj@assays$RNA@scale.data <- ges.mat 
    return(seurat.obj)
  } else {
    return(ges.mat)
  }
}

#' Performs a rank transformation on a given matrix, typically as an alternative GES generation technique.
#' If a seurat object, stores results in scale.data field; otherwise, returns rank transformation  matrix.
#'
#' @param data.object Either a Seurat object or a matrix of raw count data (genes X samples).
#' @return Rank transformed matrix, or appropriately adjusted seurat object.
#' @export
RankTransform <- function(data.object) {
  # check if seurat object
  if (class(data.object)[1] == "Seurat") {
    dat.mat <- as.matrix(seurat.obj@assays$RNA@data)
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
    seurat.obj@assays$RNA@scale.data <- rank.mat 
    return(seurat.obj)
  } else {
    return(rank.mat)
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

#' Generates a distance matrix using sqrt(1-cor(x)) as a distance metric.
#' If a seurat object, stores results in cor.dist (misc); otherwise, returns distance matrix.
#' 
#' @param data.object Either a Seurat object or a matrix of data (features X samples).
#' @param use.scaled If specified AND if data.object is a Seurat object, will use the "scale.data" field instead of "data"
#' @param cor.method Method argument for cor function; spearman by default.
#' @return Rank transformed matrix, or appropriately adjusted seurat object.
#' @export
CorDist <- function(data.object, use.scaled = FALSE, cor.method = 'spearman') {
  # check if seurat object
  if (class(data.object)[1] == "Seurat") {
    if (use.scaled) {
      dat.mat <- as.matrix(seurat.obj@assays$RNA@scale.data)
    } else {
      dat.mat <- as.matrix(seurat.obj@assays$RNA@data)
    }
  } else {
    dat.mat <- data.object
  }
  # generate distane matrix
  dist.mat <- as.dist(sqrt(1 - cor(dat.mat, method = cor.method)))
  # return
  if (class(data.object)[1] == "Seurat") {
    seurat.obj@misc[['dist.mat']] <- dist.mat 
    return(seurat.obj)
  } else {
    return(dist.mat)
  }
}
