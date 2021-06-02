#' PAM clustering over a range of k values.
#' 
#' @param dist.mat A distance matrix for the data to be clustered.
#' @param kmin The minimum k to be used (default of 2).
#' @param kmax The maximum k to be used (default of 5).
#' @param verbose Switch to print status updates.
#' @return A list of two lists; 'clusterings', which contains the cluster labels and 'sils', which has silhouette scores.
#' @export
PamKRange <- function(dist.mat, kmin = 2, kmax = 5, verbose = TRUE) { 
  ## generate clusterings for each k
  clusterings <- list()
  sils <- list()
  for (k in kmin:kmax) {
    if (verbose) { print(paste('Clustering with k=', k, '...', sep = ''))}
    # cluster
    pclust <- cluster::pam(dist.mat, k, diss = TRUE)
    clusterings[[paste('k', k, sep='')]] <- pclust
    # evaluate
    sil.score <- cluster::silhouette(pclust, dist.mat)
    sils[[paste('k', k, sep = '')]] <- mean(sil.score[,3])
  }
  return(list('clusterings' = clusterings, 'sils' = sils))
}

#' Louvain clustering over a range of resolution parameters (uses getComMembership from the MUDAN package)
#'
#' @param data.object Seurat object w/ scale.data and a distance matrix in misc$dist.mat included.
#' @param rmin Lowest resolution to try. Default of 10.
#' @param rmax Maximum resolution to try. Default of 100.
#' @param rstep Step size for resolution parameter. Default of 10.
#' @param verbose Switch to control terminal read out of progress. Default of TRUE.
#' @export
LouvainResRange <- function(data.obj, rmin = 10, rmax = 100, rstep = 10, verbose = TRUE) {
  # get matrices
  dat.mat <- data.obj@assays[[data.obj@active.assay]]@scale.data
  dist.mat <- data.obj@assays[[data.obj@active.assay]]@misc$dist.mat
  # iterate through resolution params
  res <- rmin
  clusterings <- list()
  sils <- list()
  while (res <= rmax) {
    if (verbose) {print(paste('Clustering with res=', res, '...', sep = ''))}
    # cluster
    set.seed(343)
    lclust <- MUDAN::getComMembership(t(dat.mat), k = res, method = igraph::cluster_walktrap, verbose = FALSE)
    clusterings[[paste('res', res, sep = '')]] <- lclust
    # evaluate
    sil.score <- cluster::silhouette(as.integer(lclust), dist.mat)
    sils[[paste('res', res, sep = '')]] <- mean(sil.score[,3])
    # iterate resolution param
    res <- res + rstep
  }
  # identify optimal cluster
  opt.clust <- clusterings[[which.max(sils)]]
  # add to data.object
  data.obj@assays[[data.obj@active.assay]]@misc[['pisces.cluster']] <- opt.clust
  data.obj@assays[[data.obj@active.assay]]@misc[['clustering.obj']] <- list('clusterings' = clusterings, 'sils' = sils)
  return(data.obj)
}

#' Generation of silhouette score for a list of clusters.
#'
#' @param clusterings List of clustering objects.
#' @param dist.mat Distance matrix for the data clustered in clusterings.
#' @param plotPath If specified, will save a plot of the silhouette scores.
#' @return List of silhouette scores.
#' @export
SilScoreEval <- function(clusterings, dist.mat, plotPath) {
  ## find silhouette scores for each cluster
  L <- length(clusterings)
  sil.scores <- c()
  k.vals <- c()
  for (i in 1:L) {
    # identify k for this clustering
    k <- length(table(clusterings[[i]]$clustering))
    k.vals <- c(k.vals, k)
    # find silhouette score
    sil <- cluster::silhouette(clusterings[[i]], dist.mat)
    sil.scores <- c(sil.scores, mean(sil[,3]))
  }
  ## plot if requested
  if (!missing(plotPath)) {
    plot.dat <- data.frame('k' = k.vals, 'Silhouette.Scores' = sil.scores)
    ggplot2::ggplot(plot.dat, aes(x=k, y=Silhouette.Scores)) + ggplot2::geom_point() + ggplot2::geom_line() +
      ggplot2::ggsave(plotPath, height=2, width=3)
  }
  ## return scores
  return(sil.scores)
}

#' Silhouette score optimized Louvain clustering.
#'
#' @param seurat.obj Seurat object with Neighbors detected.
#' @param dist.mat Distance matrix for use w/ Silhouette score.
#' @param rstep Size of the steps to take in resolution. Default of 0.05.
#' @param rmax Maximum resolution to try. Default of 0.5.
#' @return Seurat object with all clustering results and optimal result as active.ident.
#' @export
SSLouvain <- function(seurat.obj, dist.mat, rstep = 0.05, rmax = 0.5) {
  sil.scores <- list()
  res.val <- 0.05
  # cluster until maximum resolution exceeded
  while (res.val < rmax) {
    seurat.obj <- FindClusters(seurat.obj, resolution = res.val)
    s.clust <- as.integer(seurat.obj@active.ident) - 1
    s.score <- cluster::silhouette(s.clust, dist = dist.mat)
    sil.scores[[as.character(res.val)]] <- mean(s.score[,3])
    # increment resolution
    res.val <- res.val + rstep
  }
  # identify optimal clustering
  seurat.obj@misc[['gexp.sil.scores']] <- sil.scores
  opt.clust <- names(sil.scores)[which.max(sil.scores)]
  Idents(seurat.obj) <- seurat.obj@meta.data[[paste('SCT_snn_res.', opt.clust, sep = '')]]
  return(seurat.obj)
}