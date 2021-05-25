#' Merges two viper matrices, giving priority to one over the other.
#'
#' @param p.mat Priority viper matrix (proteins X samples). Proteins here will override those in the other matrix.
#' @param q.mat Secondary viper matrix (proteins X samples). Proteins here will fin in for gaps in the priority matrix.
#' @return A merged viper matrix.
#' @export
ViperMerge <- function(p.mat, q.mat) {
  fill.genes <- setdiff(rownames(q.mat), rownames(p.mat))
  merged.mat <- rbind(p.mat, q.mat[fill.genes,])
  return(merged.mat)
}

#' Processes ARACNe results into a regulon object compatible with VIPER.
#'
#' @param a.file ARACNe final network .tsv.
#' @param exp.mat Matrix of expression from which the network was generated (genes X samples).
#' @param out.dir Output directory for networks to be saved to.
#' @param out.name Optional argument for prefix of the file name.
#' @export
RegProcess <- function(a.file, exp.mat, out.dir, out.name = '.') {
  processed.reg <- viper::aracne2regulon(afile = a.file, eset = exp.mat, format = '3col')
  saveRDS(processed.reg, file = paste(out.dir, out.name, 'unPruned.rds', sep = ''))
  pruned.reg <- viper::pruneRegulon(processed.reg, 50, adaptive = FALSE, eliminate = TRUE)
  saveRDS(pruned.reg, file = paste(out.dir, out.name, 'pruned.rds', sep = ''))
}

#' Unwraps a nested MR list: previous functions return cluster specific master regulators as a list of lists. This funciton will unwrap that object into one, unique list.
#'
#' @param MRs List of lists, with MR names as sub-list names and MR activity as sub-list entries.
#' @param top If specified, will subset the top X regulators from each set.
#' @return Returns a de-duplicated list of MRs.
#' @export
MR_UnWrap <- function(MRs, top) {
  if (missing(top)) {
    return( unique(unlist(lapply(MRs, names), use.names = FALSE)) )
  } else {
    mr.unwrap <- lapply(MRs, function(x) {
      names(sort(x, decreasing = TRUE))[ 1:min(top, length(x)) ]
    })
    return( unique(unlist(mr.unwrap, use.names = FALSE)) )
  }
}

#' Identifies MRs for given data using stouffer integration.
#'
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param cluster Vector of cluster lables. If not included, integrates the entire matrix.
#' @param weights A named vector of sample weights. If included, stouffer integration is weighted.
#' @return Returns the stouffer integrated scores for each protien.
#' @export
StoufferMRs <- function(dat.mat, cluster, weights) {
  # generate dummy weights if missing
  if (missing(weights)) {
    weights <- as.numeric(rep(1, ncol(dat.mat))); names(weights) <- colnames(dat.mat)
  }
  # perform integration across full matrix if cluster was missing
  if (missing(cluster)) {
    sInt <- rowSums(t(t(dat.mat) * weights))
    sInt <- rowSums(t(t(dat.mat) * weights)) / sqrt(sum(weights ** 2))
    return(sort(sInt, decreasing = TRUE))
  }
  # separate cluster specific matrices
  k <- length(table(cluster))
  mrs <- list()
  for (i in 1:k) { # for each cluster
    clust.cells <- names(cluster)[which(cluster == i)]
    clust.mat <- dat.mat[, clust.cells]
    clust.weights <- weights[clust.cells]
    clust.mrs <- StoufferMRs(clust.mat, weights = clust.weights)
    mrs[[paste('c', i, sep = '')]] <- sort(clust.mrs, decreasing = TRUE)
  }
  return(mrs)
}

#' Returns the master regulators for the given data.
#'
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param method 'Stouffer' or 'ANOVA'
#' @param clustering Optional argument for a vector of cluster labels.
#' @param numMRs Number of MRs to return per cluster. Default of 50.
#' @param bottom Switch to return downregulated proteins in MR list. Default FALSE>
#' @param weights Optional argument for weights, which can be used in the Stouffer method.
#' @return Returns a list of master regulators, or a list of lists if a clustring is specified.
#' @export
GetMRs <- function(dat.mat, clustering, method, numMRs = 50, bottom = FALSE, weights, ...) {
  if (method == 'ANOVA') {
    mr.vals <- AnovaMRs(dat.mat, clustering)
  } else if (method == 'Stouffer') {
    # generate dummy weights if not specified
    if (missing(weights)) {
      weights <- rep(1, ncol(dat.mat))
      names(weights) <- colnames(dat.mat)
    }
    # recursive calls for each cluster
    if (missing(clustering)) { # no clustering specified
      mr.vals <- StoufferMRs(dat.mat, weights)
    } else {
      k <- length(table(clustering))
      mrs <- list()
      for (i in 1:k) {
        # get cluster specific matrix and weights
        clust.cells <- names(which(clustering == i))
        clust.mat <- dat.mat[, clust.cells]
        print(dim(clust.mat))
        clust.weights <- weights[clust.cells]
        # find mrs and add to list
        clust.mrs <- GetMRs(clust.mat, method = method, weights = clust.weights, numMRs = numMRs, bottom = bottom)
        print(head(clust.mrs))
        mrs[[paste('c', i, sep = '')]] <- clust.mrs
      }
      return(mrs)
    }
  } else {
    print('Invalid method: must be "Stouffer" or "ANOVA".')
  }
  # return appropriate portion of MR list
  mr.vals <- sort(mr.vals, decreasing = TRUE)
  if (bottom) {
    return(c(mr.vals[1:numMRs], tail(mr.vals, numMRs)))
  } else {
    return(mr.vals[1:numMRs])
  }
}

#' Identifies MRs on a cell-by-cell basis and returns a merged, unique list of all such MRs.
#'
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param numMRs Default number of MRs to identify in each cell. Default of 25.
#' @return Returns a list of master regulators, the unique, merged set from all cells.
#' @export
CBCMRs <- function(dat.mat, numMRs = 25) {
  # identify MRs
  cbc.mrs <- apply(dat.mat, 2, function(x) { names(sort(x, decreasing = TRUE))[1:numMRs] })
  cbc.mrs <- unique(unlist(as.list(cbc.mrs)))
  # return
  return(cbc.mrs)
}

#' Returns the viperSimilarity of the elements of two matrices.
#' 
#' @param mat.1 First matrix (features X columns).
#' @param mat.2 Second matrix (features X columns).
#' @return A matrix w/ viperSimilarity for elements of mat.1 crossed with elements of mat.2.
#' @export
TwoMatVipSim <- function(mat.1, mat.2) {
  m1.samp.num <- ncol(mat.1)
  shared.prots <- intersect(rownames(mat.1), rownames(mat.2))
  merged.mat <- cbind(mat.1[shared.prots,], mat.2[shared.prots,])
  vd.mat <- as.matrix(viper::viperSimilarity(merged.mat))
  vd.mat <- vd.mat[(m1.samp.num + 1):nrow(vd.mat), 1:m1.samp.num]
  return(vd.mat)
}