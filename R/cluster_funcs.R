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
#' @param dat.mat Matrix of data (features X samples)
#' @param dist.mat Distance matrix
#' @param rmin Lowest resolution to try. Default of 10.
#' @param rmax Maximum resolution to try. Default of 100.
#' @param rstep Step size for resolution parameter. Default of 10.
#' @param verbose Switch to control terminal read out of progress. Default of TRUE.
#' @return A list of two lists; 'clusterings', which contains the cluster labels and 'sils', which has silhouette scores
#' @export
LouvainResRange <- function(dat.mat, dist.mat, rmin = 10, rmax = 100, rstep = 10, verbose = TRUE) {
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
  return(list('clusterings' = clusterings, 'sils' = sils))
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