#' Identifies cluster specific master regulators by generating cluter-specific signatures and running NaRnEA.
#' 
#' @param norm.counts Matrix of normalized counts (features X samples).
#' @param clust.vec Clustering vector.
#' @param net.list List of ARACNe3 networks.
#' @param num.mrs Number of master regulators to return for each cluster. Default of 10.
#' @return A three member list w/ `positive` and `negative` MRs as well as the full NaRnEA results `mr.narnea`.
#' @export
cluster_signature_mrs <- function(norm.counts, clust.vec, net.list, num.mrs = 10) {
  clust.ges <- list()
  clust.names <- as.character(sort(unique(clust.vec)))
  cat('Generating cluster-specific GES...\n')
  for (cn in clust.names) {
    print(cn)
    # set sample vectors
    test.samps <- which(clust.vec == cn)
    ref.samps <- which(clust.vec != cn)
    # run wilcox test
    ges.vec <- apply(norm.counts, 1, function(x) {
      w.test <- wilcox.test(x[test.samps], x[ref.samps], alternative = "two.sided")
      rbs.cor <- 2 * w.test$statistic / (length(test.samps) * length(ref.samps)) - 1
      return(qnorm(1 - w.test$p.val) * sign(rbs.cor))
    })
    # catch nans / Inf values
    ges.vec[which(is.nan(ges.vec))] <- 0
    inf.max <- max(abs(ges.vec[which(!is.infinite(ges.vec))]))
    ges.vec[which(ges.vec == Inf)] <- inf.max + 1
    ges.vec[which(ges.vec == -Inf)] <- (-1) * inf.max - 1
    # add to list
    clust.ges[[cn]] <- ges.vec
  }
  clust.ges <- Reduce(cbind, clust.ges)
  colnames(clust.ges) <- clust.names
  # run NaRnEA
  cat('Running NaRnEA...\n')
  cluster.narnea <- meta_narnea(clust.ges, net.list)
  # pull out positive MRs
  pos.mrs <- apply(cluster.narnea$PES, 2, function(x) {
    names(sort(x, decreasing = TRUE)[1:num.mrs])
  })
  neg.mrs <- apply(cluster.narnea$PES, 2, function(x) {
    names(sort(x, decreasing = FALSE)[1:num.mrs])
  })
  # return mr object
  mr.list <- list('positive' = pos.mrs, 'negative' = neg.mrs, 'mr.narnea' = cluster.narnea)
  return(mr.list)
}


#' Identifies cluster specific master regulators using the Mann-Whitney U-test.
#' Approximates p-vals using a normal distribution for n > 30.
#' 
#' @param dat.object Seurat object w/ PISCES assay and viper matrix in 'scale.data'.
#' @param clust.vect Vector of cluster labels; if not specified, will use 'PISCES;misc$pisces.cluster' as the cluster labels
#' @return List of log p-values for pos/neg MRs in each cluster; stored in PISCES assay under 'misc' as 'mwuMRs'
#' @export
MWUMrs <- function(dat.object, clust.vect) {
  # check if seurat object
  if (class(dat.object)[1] == "Seurat") {
    dat.mat <- dat.object@assays$PISCES@scale.data
  } else {
    dat.mat <- dat.object
  }
  # get cluster vector if not specified
  if (missing(clust.vect)) {
    clust.vect <- dat.object@assays$PISCES@misc$pisces.cluster
  }
  # identify cluster names
  clust.names <- sort(unique(clust.vect))
  clust.mrs <- list()
  # cluster specific mrs
  for (cn in clust.names) {
    print(paste('Identifying MRs for cluster', cn))
    # set up labels and test statistics
    clust.samps <- names(clust.vect)[which(clust.vect == cn)]; n.1 <- length(clust.samps)
    ref.samps <- setdiff(names(clust.vect), clust.samps); n.2 <- length(ref.samps)
    u.mean <- (n.1 * n.2) / 2; u.sd <- sqrt((n.1 * n.2 * (n.1 + n.2 + 1)) / 12)
    if (n.1 < 30 | n.2 < 30) { print('WARNING: Group size <30, normal approximation may not be appropriate...') }
    # generate tests; scale; transform to p-value
    clust.wStat <- apply(dat.mat, 1, function(x) {wilcox.test(x[clust.samps], x[ref.samps])$statistic} )
    clust.zScore <- sapply(clust.wStat, function(x) {(x - u.mean) / u.sd} )
    clust.logp <- sapply(clust.zScore, function(x) {log(2) + pnorm(abs(x), lower.tail = FALSE, log.p = TRUE)})
    # check medians
    median.dif <- apply(dat.mat, 1, function(x) {sign(median(x[clust.samps]) - median(x[ref.samps]))} ) 
    # sort and return
    mr.lists <- list('positive' = sort(clust.logp[which(median.dif == 1)]),
                     'negative' = sort(clust.logp[which(median.dif == -1)]))
    clust.mrs[[as.character(cn)]] <- mr.lists
  }
  clust.mrs <- clust.mrs[sort(names(clust.mrs))]
  # add to object and return
  if (class(dat.object)[1] == "Seurat") {
    dat.object@assays$PISCES@misc[['mwuMRs']] <- clust.mrs
    return(dat.object)
  } else {
    return(clust.mrs)
  }
}

