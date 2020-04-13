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

#' Genreates breaks for a color scale based on quantiles.
#'
#' @param dat.mat Data matrix (features X samples).
#' @param n Number of breaks to generate. Default of 10.
#' @return Numeric vector of break values.
#' @export
QuantileBreaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(from = 0, to = 1, length.out = n))
  return(unique(breaks))
}

#' Creates a grid of plots for the given markers that exist in the given matrix.
#'
#' @param umap UMAP object generated from the 'UMAP' package
#' @param cluster Vector of cluster labels for each sample
#' @param pAct Matrix of protein activity (features X samples).
#' @param markers List of markers of interest.
#' @param plotTitle Title for the plot.
#' @return NULL
#' @export
MarkerGrid <- function(umap, cluster, pAct, markers, plotTitle) {
  # create plot data frame
  plot.dat <- data.frame('UMAP1' = umap$layout[,1], 'UMAP2' = umap$layout[,2],
                         'clust' = cluster)
  marker.mat <- t(pAct[ intersect(markers, rownames(pAct)) ,])
  plot.dat <- cbind(plot.dat, as.data.frame(marker.mat))
  # create plots
  clust.plot <- ggplot2::ggplot(plot.dat, ggplot2::aes(x = UMAP1, y = UMAP2, color = clust)) + ggplot2::geom_point() +
    ggplot2::ggtitle('Clusters')
  plot.list <- list('Clusters' = clust.plot)
  for (i in 4:ncol(plot.dat)) { # creat a plot for each marker and add to list
    m <- colnames(plot.dat)[i]
    m.plot <- ggplot2::ggplot(plot.dat, ggplot2::aes_string(x = 'UMAP1', y = 'UMAP2', color = m)) + ggplot2::geom_point() +
      ggplot2::ggtitle(m) + ggplot2::scale_colour_gradientn(colours = c('blue', 'white', 'red'))
    plot.list[[m]] <- m.plot
  }
  # determine dimensions
  nCol <- min(3, length(plot.list))
  nRow <- ceiling(length(plot.list) / nCol)
  # arrange plots
  marker.plot <- ggpubr::ggarrange(plotlist = plot.list, ncol = nCol, nrow = nRow)
  print(ggpubr::annotate_figure(marker.plot, top = text_grob(plotTitle, size = 24)))
}
