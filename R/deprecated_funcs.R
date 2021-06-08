#' Identifies MRs based on ANOVA analysis for a given clustering.
#'
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param clustering Clustering vector
#' @return A named vector of p-values for each protein
AnovaMRs <- function(dat.mat, clustering) {
  pVals <- c()
  group.vec <- clustering[colnames(dat.mat)]
  # perform an anova for each protein, storing pValues in a vector
  for (i in 1:nrow(dat.mat)) {
    aov.df <- data.frame('weights' = dat.mat[i,], 'group' = group.vec)
    #print(aov.df)
    aov.test <- aov(weights ~ group, aov.df)
    pVal <- summary(aov.test)[[1]][1,5]
    pVals <- c(pVals, pVal)
  }
  # name and return the vector
  names(pVals) <- rownames(dat.mat)
  return(pVals)
}

#' Performs a bootstrap t-test between two sample vectors x and y. Returns a log p-value.
#'
#' @param x Vector of test values.
#' @param y Vector of reference values.
#' @param bootstrap.num Number of bootstraps to use. Default of 100.
#' @return A signed log p-value.
LogBootstrapTTest <- function(x, y, bootstrap.num = 100) {
  x.n <- length(x); y.n <- length(y)
  log.pValue <- c()
  ## perform test for each bootstrap
  for (i in 1:bootstrap.num) {
    # create bootstraps
    x.boot <- sample(1:x.n, size = x.n, replace = TRUE)
    x.boot <- x[x.boot]
    y.boot <- sample(1:y.n, size = y.n, replace = TRUE)
    y.boot <- y[y.boot]
    # perform t.test
    test.res <- t.test(x = x.boot, y = y.boot, alternative = "two.sided")
    # generate log p-value
    log.p <- 2*pt(q = abs(test.res$statistic), df = floor(test.res$parameter), log.p = TRUE, lower.tail = FALSE)*(-sign(test.res$statistic))
    log.pValue <- c(log.pValue, log.p)
  }
  # return mean log p-value
  return(mean(log.pValue))
}

#' Performs a t-test between two sample vectors x and y. Returns a log p-value.
#' @param x Vector of test values.
#' @param y Vector of reference values.
#' @param bootstrap.num Number of bootstraps to use. Default of 100.
#' @return A signed log p-value.
LogTTest <- function(x, y) {
  test.res <- t.test(x, y, alternative = 'two.sided')
  log.p <- 2*pt(q = abs(test.res$statistic), df = floor(test.res$parameter), log.p = TRUE, lower.tail = FALSE)*(-sign(test.res$statistic))
  return(log.p)
}

#' Identifies MRs based on a bootstraped Ttest between clusters.
#'
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param class.vect Vector of cluster labels / class vectors.
#' @param bootstrap.num Number of bootstraps to use. Default of 10
#' @return Returns a list of lists; each list is a vector of sorted log p-values for each cluster.
BTTestMRs <- function(dat.mat, class.vect, bootstrap.num = 100) {
  # get class names
  class.names <- unique(class.vect)
  mrs <- list()
  # iterate thrugh classes
  for (cn in class.names) { 
    print(paste('Identifying MRs for class ', cn, '...', sep = ''))
    mrs.mat <- matrix(0L, nrow = nrow(dat.mat), ncol = bootstrap.num)
    rownames(mrs.mat) <- rownames(dat.mat)
    # split test and ref matrices
    cn.vect <- which(class.vect == cn)
    test.mat <- dat.mat[, cn.vect]; ref.mat <- dat.mat[, -cn.vect]
    t.n <- ncol(test.mat); r.n <- ncol(ref.mat)
    # for each bootstrap
    for (b in 1:bootstrap.num) {
      test.boot <- test.mat[, sample(colnames(test.mat), size = t.n, replace = TRUE)]
      ref.boot <- ref.mat[, sample(colnames(ref.mat), size = t.n, replace = TRUE)]
      # for each gene
      for (g in rownames(dat.mat)) {
        mrs.mat[g, b] <- LogTTest(test.boot[g,], ref.boot[g,])
      }
    }
    # sort and add to list
    mList <- sort(rowMeans(mrs.mat), decreasing = TRUE)
    mrs[[cn]] <- mList
  }
  # return
  return(mrs)
}

#' Filters out genes with no expression and low quality cells.
#'
#' @param raw.mat Matrix of raw gene expression data (genes X samples).
#' @param minCount Minimum number of reads in a cell. Default of 1000.
#' @param maxCount Maximum number of reads in a cell. Default of 100000.
#' @param minGeneReads Minimum number of reads for a gene to be kept. Default of 1 (any gene with no reads will be removed).
#' @return Quality controlled matrix.
QCTransform <- function(raw.mat, minCount = 1000, maxCount = 100000, minGeneReads = 1) {
  filt.mat <- raw.mat[, colSums(raw.mat) > minCount & colSums(raw.mat) < maxCount]
  filt.mat <- filt.mat[ rowSums(filt.mat) >= minGeneReads ,]
  rem.genes <- nrow(raw.mat) - nrow(filt.mat); rem.cells <- ncol(raw.mat) - ncol(filt.mat)
  print(paste('Removed ', rem.genes, ' gene(s) and ', rem.cells, ' cell(s).', sep =''))
  return(filt.mat)
}

#' Returns vector of mitochondrial percentages for the given samples.
#'
#' @param dat.mat Matrix of raw gene expression data (genes X samples).
#' @param species
#' @return Returns named vector of mitochondrial gene percentage.
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

#' Identifies MRs for given data using stouffer integration.
#'
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param cluster Vector of cluster lables. If not included, integrates the entire matrix.
#' @param weights A named vector of sample weights. If included, stouffer integration is weighted.
#' @return Returns the stouffer integrated scores for each protien.
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