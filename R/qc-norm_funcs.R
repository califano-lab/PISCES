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

#' Returns vector of mitochondrial percentages for the given samples.
#'
#' @param dat.mat Matrix of raw gene expression data (genes X samples).
#' @param species
#' @retun Returns named vector of mitochondrial gene percentage.
#' @export
MTPercent <- function(dat.mat, mt.genes) {
  mt.count <- colSums(dat.mat[ intersect(rownames(dat.mat), mt.genes) ,])
  total.count <- colSums(dat.mat)
  mt.percent <- mt.count / total.count
  head(mt.percent)
  return( mt.percent )
}

#' Filters data based on percentage of mitochondrial gens.
#'
#' @param raw.mat Matrix of raw gene expression data (genes X samples)
#' @param mt.genes Path to .csv file with ENSG and Hugo names for mitochondrial genes.
#' @param mt.thresh Threshold above which cells will be removed. Default of 0.15
#' @export
MTFilter <- function(dat.mat, mt.genes, mt.thresh = 0.1) {
  ## find mt percentages
  mt.perc <- MTPercent(dat.mat, mt.genes)
  ## filter matrix
  thresh.cells <- names(mt.perc)[which(mt.perc < mt.thresh)]
  rem.cells <- ncol(dat.mat) - length(thresh.cells)
  print(paste('Removed', rem.cells, 'cell(s) with too many MT reads', sep = ' '))
  dat.mat <- dat.mat[, thresh.cells ]
  return(dat.mat)
}

#' Filters out genes with no expression and low quality cells.
#'
#' @param raw.mat Matrix of raw gene expression data (genes X samples).
#' @param minCount Minimum number of reads in a cell. Default of 1000.
#' @param maxCount Maximum number of reads in a cell. Default of 100000.
#' @param minGeneReads Minimum number of reads for a gene to be kept. Default of 1 (any gene with no reads will be removed).
#' @return Quality controlled matrix.
#' @export
QCTransform <- function(raw.mat, minCount = 1000, maxCount = 100000, minGeneReads = 1) {
  filt.mat <- raw.mat[, colSums(raw.mat) > minCount & colSums(raw.mat) < maxCount]
  filt.mat <- filt.mat[ rowSums(filt.mat) >= minGeneReads ,]
  rem.genes <- nrow(raw.mat) - nrow(filt.mat); rem.cells <- ncol(raw.mat) - ncol(filt.mat)
  print(paste('Removed ', rem.genes, ' gene(s) and ', rem.cells, ' cell(s).', sep =''))
  return(filt.mat)
}

#' Performas a CPM normalization on the given data.
#'
#' @param dat.mat Matrix of gene expression data (genes X samples).
#' @param l2 Optional log2 normalization switch. Default of False.
#' @param pseudo Optional pseudo count logical. Default of False.
#' @return Returns CPM normalized matrix
#' @export
CPMTransform <- function(dat.mat, l2 = FALSE, pseudo = FALSE) {
  if (pseudo) {
    dat.mat <- dat.mat + 1
  }
  cpm.mat <- t(t(dat.mat) / (colSums(dat.mat) / 1e6))
  if (l2) {
    cpm.mat <- log2(cpm.mat + 1)
  }
  return(cpm.mat)
}

#' Performs a rank transformation on a given matrix.
#'
#' @param dat.mat Matrix of data, usually gene expression (genes X samples).
#' @return Rank transformed matrix.
#' @export
RankTransform <- function(dat.mat) {
  rank.mat <- apply(dat.mat, 2, rank)
  median <- apply(rank.mat, 1, median)
  mad <- apply(rank.mat, 1, mad)
  rank.mat <- (rank.mat - median) / mad
  return(rank.mat)
}

#' Generates a gene expression signature (GES) using internal normalization.
#'
#' @param cpm.mat Matrix of CPM-normalized gene expression data (genes X samples).
#' @return GES matrix.
#' @export
GESTransform <- function(cpm.mat) {
  ges.mat <- t(apply(cpm.mat, 1, function(x) {
    (x - mean(x)) / sd(x)
  }))
  return(ges.mat)
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
