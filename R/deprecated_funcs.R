#' Identifies MRs based on ANOVA analysis for a given clustering.
#'
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param clustering Clustering vector
#' @return A named vector of p-values for each protein
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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