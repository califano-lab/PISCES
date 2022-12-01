#' Identifies cluster specific master regulators by generating cluter-specific signatures and running NaRnEA.
#' 
#' @param norm.counts Matrix of normalized counts (features X samples).
#' @param clust.vec Clustering vector.
#' @param net.list List of ARACNe3 networks.
#' @return A list; table of positive MRs (`positive`) and negative MRs (`negative`); 
#' the cluster-specific NaRnEA results (`cluster.narnea`).
#' @export
cluster_signature_mrs <- function(norm.counts, clust.vec, net.list) {
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
  # identify significant Mrs in each cluster and sort into positive/negative sets
  sig.prots <- sig_features(cluster.narnea$NES)
  pos.mrs <- lapply(clust.names, function(x) {
    pos.sig.clust.prots <- intersect(sig.prots[[x]], names(which(cluster.narnea$PES[,x] > 0)))
    pos.sig.pes <- cluster.narnea$PES[pos.sig.clust.prots, x]
    names(pos.sig.pes) <- pos.sig.clust.prots
    return(names(sort(pos.sig.pes, decreasing = TRUE)))
  })
  neg.mrs <- lapply(clust.names, function(x) {
    neg.sig.clust.prots <- intersect(sig.prots[[x]], names(which(cluster.narnea$PES[,x] < 0)))
    neg.sig.pes <- cluster.narnea$PES[neg.sig.clust.prots, x]
    names(neg.sig.pes) <- neg.sig.clust.prots
    return(names(sort(neg.sig.pes, decreasing = TRUE)))
  })
  names(pos.mrs) <- clust.names; names(neg.mrs) <- clust.names
  # return mr object
  mr.list <- list('positive' = pos.mrs, 'negative' = neg.mrs, 'cluster.narnea' = cluster.narnea)
  return(mr.list)
}

#' Identified MRs based on a pairwise comparison of clusters.
#' 
#' @param norm.counts Matrix of normalized counts (features X samples).
#' @param clust.vec Clustering vector.
#' @param net.list List of ARACNe3 networks.
#' @param num.mrs Number of master regulators to return for each cluster. Default of 10.
#' @return A three member list w/ `positive` and `negative` MRs as well as the full NaRnEA results `mr.narnea`.
#' @export
pairwise_cluster_mrs <- function(norm.counts, clust.vec, net.list, num.mrs = 10) {
  # generate pairwise GES
  ges.vecs <- list()
  clust.names <- as.character(sort(unique(clust.vec)))
  num.clust <- length(clust.names)
  cat("Generating pairwise cluster GES...\n")
  for (i in 1:(num.clust - 1)) {
    for (j in (i+1):num.clust) {
      i.clust <- clust.names[i]
      j.clust <- clust.names[j]
      # get samples
      i.samps <- which(clust.vec == i.clust)
      j.samps <- which(clust.vec == j.clust)
      # generate GES
      ges.vec <- apply(norm.counts, 1, function(x) {
        w.test <- wilcox.test(x[i.samps], x[j.samps], alternative = "two.sided")
        rbs.cor <- 2 * w.test$statistic / (length(i.samps) * length(j.samps)) - 1
        return(qnorm(1 - w.test$p.val) * sign(rbs.cor))
      })
      # catch nans / Inf values
      ges.vec[which(is.nan(ges.vec))] <- 0
      inf.max <- max(abs(ges.vec[which(!is.infinite(ges.vec))]))
      ges.vec[which(ges.vec == Inf)] <- inf.max + 1
      ges.vec[which(ges.vec == -Inf)] <- (-1) * inf.max - 1
      # add to list
      ges.vecs[[paste(i.clust, 'v', j.clust, sep = '.')]] <- ges.vec
    }
  }
  clust.ges <- Reduce(cbind, ges.vecs)
  colnames(clust.ges) <- names(ges.vecs)
  # run NaRnEA
  cat('Running NaRnEA...\n')
  pairwise.narnea <- meta_narnea(clust.ges, net.list)
  # pull out positive MRs
  pos.mrs <- apply(pairwise.narnea$PES, 2, function(x) {
    names(sort(x, decreasing = TRUE)[1:num.mrs])
  })
  colnames(pos.mrs) <- names(ges.vecs)
  neg.mrs <- apply(pairwise.narnea$PES, 2, function(x) {
    names(sort(x, decreasing = FALSE)[1:num.mrs])
  })
  colnames(neg.mrs) <- names(ges.vecs)
  # return mr object
  mr.list <- list('positive' = pos.mrs, 'negative' = neg.mrs, 'cluster.narnea' = cluster.narnea)
  return(mr.list)  
}

#' Returns an MR object with MRs selected via Stouffer integration.
#' 
#' @param narnea.res NaRnEA result object with `PES` and `NES` members.
#' @param clust.vec Vector of categorical labels.
#' @param p.thresh p-value threshold for calling significant master regulators. Default of 0.05.
#' @return A list; table of positive MRs (`positive`) and negative MRs (`negative`); 
#' Stouffer integrated NaRnEA results by cluster (`cluster.narnea`).
#' @export
stouffer_cluster_mrs <- function(narnea.res, clust.vec, p.thresh = 0.05) {
  # perform stoufferintegration
  cluster.stouffer <- stouffer_narnea(narnea.res, clust.vec)
  clust.names <- colnames(cluster.stouffer$PES)
  # identify significant Mrs in each cluster and sort into positive/negative sets
  sig.prots <- sig_features(cluster.stouffer$NES)
  pos.mrs <- lapply(clust.names, function(x) {
    pos.sig.clust.prots <- intersect(sig.prots[[x]], names(which(cluster.stouffer$PES[,x] > 0)))
    pos.sig.pes <- cluster.stouffer$PES[pos.sig.clust.prots, x]
    names(pos.sig.pes) <- pos.sig.clust.prots
    return(names(sort(pos.sig.pes, decreasing = TRUE)))
  })
  neg.mrs <- lapply(clust.names, function(x) {
    neg.sig.clust.prots <- intersect(sig.prots[[x]], names(which(cluster.stouffer$PES[,x] < 0)))
    neg.sig.pes <- cluster.stouffer$PES[neg.sig.clust.prots, x]
    names(neg.sig.pes) <- neg.sig.clust.prots
    return(names(sort(neg.sig.pes, decreasing = TRUE)))
  })
  names(pos.mrs) <- clust.names; names(neg.mrs) <- clust.names
  # return in mr list format
  return(list('positive' = pos.mrs, 'negative' = neg.mrs, 'cluster.narnea' = cluster.stouffer))
}

#' Identifies cluster master regulators based on Cohen's Kappa.
#' 
#' @param narnea.res NaRnEA result object with `PES` and `NES` members.
#' @param clust.vec Vector of clustering labels.
#' @param kappa.thresh Threshold to consider a protein an MR. Default of 0.5.
#' @return A list; table of positive MRs (`positive`) and negative MRs (`negative`); 
#' Stouffer integrated NaRnEA results by cluster (`cluster.narnea`); Kappa values for each protein in each cluster (`kappa.df`).
#' @export
kappa_cluster_mrs <- function(narnea.res, clust.vec, kappa.thresh = 0.5) {
  # perform stouffer integration
  cluster.stouffer <- stouffer_narnea(narnea.res, clust.vec)
  clust.names <- colnames(cluster.stouffer$PES)
  # process cluster information, then find kappa values for each cluster
  clust.table <- table(clust.vec)
  c.vec.base <- rep('not.clust', length(clust.vec)); names(c.vec.base) <- names(clust.vec)
  kappa.lists <- list()
  for (cn in names(clust.table)) {
    cat(paste("Calculating MRs for ", cn, "...\n", sep = ''))
    # construct vector of this cluster versus all
    c.vec <- c.vec.base
    c.vec[which(clust.vec == cn)] <- cn
    # subsample if too long; [TODO: more robust sampling]
    if (length(c.vec) > 500) {c.vec <- c.vec[sample(names(c.vec), 500)]}
    # calculate cohen kappa for each feature
    feature.kappa <- apply(narnea.res$PES[, names(c.vec)], 1, function(x) {
      opt_kappa(c.vec, x)
    })
    # calculate mean pes
    kappa.lists[[as.character(cn)]] <- feature.kappa
  }
  kappa.df <- Reduce(cbind, kappa.lists)
  colnames(kappa.df) <- clust.names
  # create pos / neg mr tables
  pos.mrs <- lapply(clust.names, function(x) {
    pos.prots <- names(which(cluster.stouffer$PES[,x] > 0))
    kappa.sort <- sort(kappa.df[,x], decreasing = TRUE)
    kappa.sort <- kappa.sort[which(kappa.sort > kappa.thresh)]
    pos.prots.intersect <- intersect(names(kappa.sort), pos.prots)
    return(pos.prots.intersect)
  })
  neg.mrs <- lapply(clust.names, function(x) {
    neg.prots <- names(which(cluster.stouffer$PES[,x] < 0))
    kappa.sort <- sort(kappa.df[,x], decreasing = TRUE)
    kappa.sort <- kappa.sort[which(kappa.sort > kappa.thresh)]
    neg.prots.intersect <- intersect(names(kappa.sort), neg.prots)
    return(neg.prots.intersect)
  })
  names(pos.mrs) <- clust.names
  names(neg.mrs) <- clust.names
  # return mr object
  return(list('positive' = pos.mrs, 'negative' = neg.mrs, 'cluster.narnea' = cluster.stouffer,
              'kappa.df' = kappa.df))
}

#' Identifies cluster-specific master regulators using a Kruskal-Wallis test.
#' 
#' @param narnea.res NaRnEA result object with `PES` and `NES` members.
#' @param clust.vec Vector of categorical labels.
#' @param p.thresh p-value threshold for calling significant master regulators. Default of 0.05.
#' @return A list; table of positive MRs (`positive`) and negative MRs (`negative`); 
#' Kruskal-Wallis results (`kw.test`) with log p-values and group residuals..
#' @export
kw_cluster_mrs <- function(narnea.res, clust.vec, p.thresh = 0.05) {
  # run kw on each protein
  kw.test <- apply(narnea.res$PES, 1, function(x) {
    log_kw(x, clust.vec)
  })
  log.p.vec <- sapply(kw.test, function(x) {x$log.p})
  residual.mat <- Reduce('rbind', lapply(kw.test, function(x) {x$group.residuals}))
  rownames(residual.mat) <- names(kw.test)
  # determine significant proteins
  adj.p <- p.adjust(exp(log.p.vec), method = 'BH')
  sig.prots <- names(which(adj.p < p.thresh))
  # assign significant proteins to clusters based on residuals
  pos.mrs <- apply(residual.mat[sig.prots,], 2, function(x) {
    return(names(sort(x[which(x > 0)], decreasing = TRUE)))
  })
  neg.mrs <- apply(residual.mat[sig.prots,], 2, function(x) {
    return(names(sort(x[which(x < 0)], decreasing = FALSE)))
  })
  # return 
  return(list('positive' = pos.mrs, 'negative' = neg.mrs, 'kw.test' = kw.test))
}

#' Determines the optimal possible Cohen's Kappa that can be generated by the given feature.
#' 
#' @param class.vec Vector of categorical labels.
#' @param feature.vec Vector of numerical feature values. Assumes the same ordering as `class.vec`.
#' @return THe optimal possible value of Cohen's Kappa for a classifier generated using this feature.
#' @export
opt_kappa <- function(class.vec, feature.vec) {
  # process into data frame sorted by feature vector
  n.samps <- length(class.vec)
  dat.df <- data.frame('class' = class.vec, 'feature' = feature.vec)
  dat.df <- dat.df[order(dat.df$feature),]
  # find candidate break points
  break.points <- which(dat.df$class[2:n.samps] != dat.df$class[1:(n.samps - 1)]) + 1
  break.vals <- dat.df$feature[break.points]
  # find the optimal kappa; break point itself is not relevant
  opt.kappa <- -Inf
  for (bv in break.vals) {
    break.class <- (dat.df$feature < bv)
    b.kappa <- cohen_kappa(break.class, dat.df$class)
    if (b.kappa > opt.kappa) { opt.kappa <- b.kappa }
  }
  # return optimal kappa
  return(opt.kappa)
}

#' Calculates the Cohen's Kappa between two categorical vectors.
#' 
#' @param vec.a First vector of categorical labels.
#' @param vec.b Second vector of categorical labels. Assumes same ordering as `vec.a`.
#' @return Cohen's Kappa
#' @export
cohen_kappa <- function(vec.a, vec.b) {
  # build confusion matrix; account for single-group classification vector
  conf.table <- table(vec.a, vec.b)
  if (nrow(conf.table) == 1) {conf.table <- rbind(conf.table, c(0,0))}
  tot.samps <- length(vec.a)
  # calculate kappa
  agree.val <- sum(diag(conf.table)) / tot.samps
  agree.chance <- (sum(conf.table[,1]) * sum(conf.table[1,]) + 
                     sum(conf.table[,2]) * sum(conf.table[2,])) / tot.samps**2
  k.val <- (agree.val - agree.chance) / (1 - agree.chance)
  # return
  return(k.val)
}

#' Returns a list of significant proteins for each sample (column) in the given NES matrix.
#' Identifies significant proteins by Benjamini-Hochberg at the given significance level.
#' 
#' @param nes.mat Matrix of NES scores (features X columns), such as that generated by NaRnEA analysis.
#' @param p.val.thresh p-value threshold. Default of 0.05.
#' @return List of lists; one list for each sample (column), containing the names of significant proteins.
#' @export
sig_features <- function(nes.mat, p.val.thresh = 0.05) {
  sig.feats <- apply(nes.mat, 2, function(x) {
    p.vals <- p.adjust(pnorm(abs(x), lower.tail = FALSE), method = 'BH')
    return(names(which(p.vals < p.val.thresh)))
  })
  names(sig.feats) <- colnames(nes.mat)
  return(sig.feats)
}

#' Returns a log p-value from a Kruskal-Wallis test for the given data.
#' 
#' @param feat.vec Vector of values to be tested.
#' @param group.vec Vector of group labels. Assumed to be in the same order as `feat.vec`.
#' @return A list; `log.p` containing the log p-value and `group.residuals` (group mean rank - overall mean rank).
#' @export
log_kw <- function(feat.vec, group.vec) {

  # calculate overall values
  rank.feat <- rank(feat.vec)
  group.table <- table(group.vec)
  group.names <- names(group.table)
  num.groups <- length(group.names)
  num.samps <- length(feat.vec)
  mean.rank <- 0.5 * (num.samps + 1)
  
  # calculate group statistics
  group.mean.rank <- sapply(group.names, function(x) {
    return(mean(rank.feat[which(group.vec == x)]))
  })
  names(group.mean.rank) <- group.names
  
  # calculate h-statistic
  h.num <- sum(group.table * (group.mean.rank - mean.rank)**2)
  h.dem <- sum((rank.feat - mean.rank)**2)
  h.val <- (num.samps - 1) * h.num / h.dem
  
  # calculate p-value
  log.p <- pchisq(h.val, df = num.groups - 1, lower.tail = FALSE, log.p = TRUE)
  
  # return
  return(list('log.p' = log.p, 'group.residuals' = (group.mean.rank - mean.rank)))
}

#' Identifies cluster specific master regulators using the Mann-Whitney U-test.
#' Approximates p-vals using a normal distribution for n > 30.
#' 
#' @param dat.object Seurat object w/ PISCES assay and viper matrix in 'scale.data'.
#' @param clust.vect Vector of cluster labels; if not specified, will use 'PISCES;misc$pisces.cluster' as the cluster labels
#' @return List of log p-values for pos/neg MRs in each cluster; stored in PISCES assay under 'misc' as 'mwuMRs'
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

#' Returns set of candidate master regulators or markers from the given mr.list object.
#' 
#' @param mr.list MR List, such as that generated by `stouffer_cluster_mrs` or `kappa_cluster_mrs`.
#' @param num.mrs Number of MRs to return for each cluster. Defualt of 10.
#' @param reg.class Class of proteins to return. One of `c('regulator', 'marker')`, with `regulator` as default.
#' @param reg.sign Sign of proteins to return. One of `c('pos', 'neg')`, with `pos` as default.
#' @return A dataframe with proteins and assocaited clusters.
#' @export
get_mr_set <- function(mr.list, num.mrs = 10, reg.class = c('regulator', 'marker'), reg.sign = c('pos', 'neg')) {
  require(stringr)
  
  match.arg(reg.class)
  if (missing(reg.class)) {reg.class = 'regulator'}
  match.arg(reg.sign)
  if (missing(reg.sign)) {reg.sign = 'pos'}
  
  # detect no regulators presnet
  mr.unlist <- switch(reg.sign,
                      'pos' = unlist(mr.list$positive),
                      'neg' = unlist(mr.list$negative))
  if (length(mr.unlist) == 0) {
    cat("Error: No regulators of the selected sign in the given list, possibly because none are significant.\n")
    return(NULL)
  }
  
  # detect species and load regulator set
  if (str_detect(unlist(mr.list$positive)[1], "[[:lower:]]")) { # mouse 
    data("mouse_regulators")
    reg.list <- mouse_regulators
  } else { # human
    data("human_regulators")
    reg.list <- human_regulators
  }
  
  # set reg set
  if (reg.class == 'regulator') {
    reg.inds <- which(names(reg.list) != 'surf')
  } else {
    reg.inds <- which(names(reg.list) == 'surf')
  }
  
  # select proteins
  if (reg.sign == 'pos') {
    clust.regs <- mr.list$positive
  } else {
    clust.regs <- mr.list$negative
  }
  
  # generate protein df and return
  reg.set <- unlist(reg.list[reg.inds])
  prot.lists <- lapply(clust.regs, function(x) {
    int.regs <- intersect(x, reg.set)
    if (length(int.regs) == 0) { return(NULL) }
    return(int.regs[1:min(num.mrs, length(int.regs))])
  })
  prot.df <- data.frame('protein' = unlist(prot.lists),
                        'cluster' = rep(names(prot.lists), sapply(prot.lists, length)))
  
  return(prot.df)
}

