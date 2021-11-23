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

#' Runs Louvain clusering on the given object with the specified KNN graph.
#' 
#' @param data.obj Seurat object w/ dist.mat object in misc of active assay.
#' @param num.neighbors Number of neighbors to use in the KNN graph. Default of 5.
#' @return Vector of cluster assignments.
#' @export
LouvainClust <- function(data.obj, num.neighbors = 5) {
  # get KNN matrix
  knn.mat <- KNN(data.obj@assays[[data.obj@active.assay]]@misc$dist.mat,
                 num.neighbors)
  # generate graph
  adj.mat <- matrix(0L, nrow = nrow(knn.mat), ncol = nrow(knn.mat))
  colnames(adj.mat) <- rownames(knn.mat)
  rownames(adj.mat) <- rownames(knn.mat)
  for (i in 1:nrow(knn.mat)) {
    adj.mat[i, knn.mat[i,]] <- 1
  }
  graph.obj <- graph_from_adjacency_matrix(adj.mat)
  graph.obj <- as.undirected(graph.obj)
  # create clustering
  l.clust <- unlist(as.list(membership(cluster_louvain(graph.obj))))
  return(l.clust)
}

#' Runs Louvain clustering w/ each of the specified number of neighbors used in the KNN graph.
#' 
#' @param data.obj Seurat object w/ dist.mat object in misc. of active assay.
#' @param kmin Minimum number of neighbors to use. Default of 5.
#' @param kmax Maximum number of neighbors to use. Default of 50.
#' @param kstep Steps to take between values of k. Default of 5.
#' @return Updated data.obj w/ clustering.obj and pisces.cluster added to active assay.
#' @export
LouvainKRange <- function(data.obj, kmin = 5, kmax = 50, kstep = 5) {
  # create lists
  cluster.list <- list()
  sil.list <- list()
  # grab distance matrix
  dist.mat <- data.obj@assays[[data.obj@active.assay]]@misc$dist.mat
  # setup iteration
  k <- kmin
  while (k <= kmax) {
    print(paste("Clustering with k = ", k, "...", sep = ''))
    # generate clustering and silhouette score
    clust.vec <- LouvainClust(data.obj, k)
    sil.score <- cluster::silhouette(clust.vec, dist.mat)
    # add to list
    k.ind <- paste('k', k, sep = '.')
    cluster.list[[k.ind]] <- clust.vec
    sil.list[[k.ind]] <- mean(sil.score[,3])
    # iterate
    k <- k + kstep
  }
  # identify optimal cluster
  opt.clust <- cluster.list[[which.max(sil.list)]]
  # add to data.object
  data.obj@assays[[data.obj@active.assay]]@misc[['pisces.cluster']] <- opt.clust
  data.obj@assays[[data.obj@active.assay]]@misc[['clustering.obj']] <- list('clusterings' = cluster.list, 'sils' = sil.list)
  return(data.obj)
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
    print(head(lclust))
    # evaluate
    if (length(unique(as.integer(lclust)) == 0)) {
      print("singular clustering")
      sils[[paste('res', res, sep = '')]] <- NA
    } else {
      sil.score <- cluster::silhouette(as.integer(lclust), dist.mat)
      sils[[paste('res', res, sep = '')]] <- mean(sil.score[,3])
    }
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