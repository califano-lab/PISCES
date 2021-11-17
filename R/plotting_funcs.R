#' Identifies a vector of colors for a given number of clusters. Internal function.
#'
#' @param k Number of clusters.
#' @param offset Optional argument to shift colors along color wheel.
#' @return A vector of hues
#' @export
ClusterColors <- function(k, offset = 0) {
  hues <- seq(15, 375, length = k + 1) + offset
  return(hcl(h = hues, l = 65, c = 100)[1:k])
}

#' Generates breaks for a color scale based on quantiles.
#'
#' @param dat.mat Data matrix (features X samples).
#' @param n Number of breaks to generate. If not specified, uses first three stdevs.
#' @return Numeric vector of break values.
#' @export
QuantileBreaks <- function(xs, n) {
  if (!missing(n)) {
    breaks <- quantile(xs, probs = seq(from = 0, to = 1, length.out = n))
  } else {
    breaks <- quantile(xs, c(0.003, 0.05, 0.32, 0.5, 0.68, 0.95, 0.997))
  }
  return(unique(breaks))
}

#' Creates a grid of plots for the given markers that exist in the given matrix.
#'
#' @param data.object Seurat object w/ PISCES assay, VIPER, and misc$umap
#' @param marker.list List of gene names w/ markers of interest.
#' @return Plot object from ggarrange
#' @export
MarkerGrid <- function(vip.obj, marker.list) {
  # prune marker.list
  marker.list <- intersect(marker.list, rownames(vip.obj@assays$PISCES@scale.data))
  # make plot.dat
  plot.dat <- data.frame('UMAP1' = vip.obj@assays$PISCES@misc$umap[,1],
                         'UMAP2' = vip.obj@assays$PISCES@misc$umap[,2])
  plot.dat <- cbind(plot.dat, t(vip.obj@assays$PISCES@scale.data[marker.list,]))
  # make plots
  plot.list <- list()
  for (p in colnames(plot.dat)[3:ncol(plot.dat)]) {
    plot.list[[p]] <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes_string(color = p)) +
      ggtitle(p) + scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')
  }
  # arrange
  nCol <- min(3, length(plot.list))
  nRow <- ceiling(length(plot.list) / nCol)
  # plot object
  plot.obj <- ggpubr::ggarrange(plotlist = plot.list, ncol = nCol, nrow = nRow)
  return(plot.obj)
}

#' Returns color gradient for the specified data type (green/purple for Gene Expression; red/blue for VIPER)
#'
#' @param num.breaks Number of breaks in the gradient.
#' @param data.type Type of data to use; either 'gexp' or 'vip'
#' @return Vector of colors.
#' @export
ColorLevels <- function(num.colors, data.type) {
  if (data.type == 'gexp') {
    col.func <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'PRGn')))
  } else if (data.type == 'vip') {
    col.func <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'RdBu')))
  } else {
    print("Error: Not a valid data type; must be one of 'gexp' or 'vip'")
    return(0)
  }
  return(col.func(num.colors))
}

#' Generates basic quality control plots from raw gene expression data.
#'
#' @param seurat.obj A Seurat object w/ MT% added to metadata.
#' @export
QCPlots <- function(seurat.obj) {
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
  ## make final plot
  plot.obj <- ggpubr::ggarrange(plotlist = list(p1, p2, p3), ncol = 3)
  return(plot.obj)
}

#' Generates a heatmap of the cluster-specific master regulators.
#' 
#' @param pisces.obj Seurat object w/ PISCES assay, VIPER matrix, and master regulators.
#' @param num.mrs Number of master regulators to include for each plot. Default of 10.
#' @param clust.vect Optional argument for cluster labels; if not included, uses "PISCES;misc$pisces.cluster"
#' @param plot.title Optional plot title argument. If not specified, no title included.
#' @return Returns the pheatmap plot object.
#' @export
MRHeatmap <- function(pisces.obj, num.mrs = 10, clust.vect, plot.title = '') {
  # check for PISCES assay
  if (!HasPISCESAssay(pisces.obj)) {
    print("Error: No PISCES assay in supplied object...")
    return(NULL)
  }
  # fetch objects
  vip.mat <- pisces.obj@assays$PISCES@scale.data
  vip.mrs <- pisces.obj@assays$PISCES@misc$mwuMRs
  vip.mrs <- vip.mrs[order(names(vip.mrs))]
  if (missing(clust.vect)) {
    clust.vect <- pisces.obj@assays$PISCES@misc$pisces.cluster
  }
  # build mr set
  mr.set <- unique(unlist(lapply(vip.mrs, function(x) {names(x$positive[1:num.mrs])} )))
  # build plot matrix and color breaks
  cell.order <- names(sort(clust.vect))
  plot.mat <- vip.mat[mr.set, cell.order]
  mat.breaks <- QuantileBreaks(plot.mat, 100)
  # cluster annotation
  annot.df <- data.frame('Cluster' = as.factor(clust.vect))
  clust.colors <- ClusterColors(length(unique(clust.vect))); names(clust.colors) <- sort(unique(clust.vect))
  annot.color <- list('Cluster' = clust.colors)
  print(dim(plot.mat))
  # make plot
  plot.obj <- pheatmap::pheatmap(plot.mat, main = plot.title, 
                                 annotation_col = annot.df, annotation_colors = annot.color,
                                 cluster_cols = FALSE, show_colnames = FALSE,
                                 cluster_rows = FALSE, show_rownames = TRUE,
                                 breaks = mat.breaks, color = ColorLevels(length(mat.breaks) - 1, 'vip'))
  return(plot.obj)
}

#' Generates a scatter plot of the clusters; must have 'umap' object in 'misc'
#' 
#' @param pisces.obj Seurat object w/ PISCES assay, VIPER matrix, and master regulators.
#' @param clust.vect Optional argument for cluster labels; if not included, uses "PISCES;misc$pisces.cluster"
#' @return Returns ggplot object
#' @export
PlotClusters <- function(pisces.obj, clust.vect) {
  # get clustering from object
  if (missing(clust.vect)) {
    clust.vect <- pisces.obj@assays$PISCES@misc$pisces.cluster
  }
  # make plot.dat
  plot.dat <- data.frame('UMAP1' = pisces.obj@assays$PISCES@misc$umap[,1],
                         'UMAP2' = pisces.obj@assays$PISCES@misc$umap[,2],
                         'Cluster' = as.factor(clust.vect))
  # return plot
  clust.plot <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cluster))
  return(clust.plot)
}
