#' Stouffer integrates the given vector of data.
#' 
#' @param dat.vect Vector of data to be integrated
#' @param weight.vect Vector of weights, if specified.
#' @return Stouffer integrated value for the given data.
#' @export
StoufferIntegrate <- function(dat.vect, weight.vect) {
  if (!missing(weight.vect)) {
    s.int <- sum(dat.vect * weight.vect) / sqrt(sum(weight.vect**2))
  } else {
    s.int <- sum(dat.vect) / length(dat.vect)
  }
  return(s.int)
}

#' Master regulators by Stouffer integration.
#' 
#' @param dat.mat Matrix of data (features X samples).
#' @param clust.vect Vector of cluster labels.
#' @param weight.vect Vector of weights, if specified.
#' @return Integrated 
StoufferMRs <- function(dat.mat, clust.vect) {
  
}

#' Identifies cluster specific master regulators using the Mann-Whitney U-test.
#' Approximates p-vals using a normal distribution for n > 30.
#' 
#' @param dat.mat Matrix of data.
#' @param clust.vect Vector of cluster labels.
#' @return For each cluster, two sorted lists of log p-values, split by positive / negative
#' @export
MWUMrs <- function(dat.mat, clust.vect) {
  # identify cluster names
  clust.names <- unique(clust.vect)
  clust.mrs <- list()
  # cluster specific mrs
  for (cn in clust.names) {
    print(paste('Identifying MRs for cluster', cn))
    # set up labels and test statistics
    clust.samps <- which(clust.vect == cn); n.1 <- length(clust.samps)
    ref.samps <- which(clust.vect != cn); n.2 <- length(ref.samps)
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
    clust.mrs[[cn]] <- mr.lists
  }
  return(clust.mrs)
}

#' Writes a given master regulator object to tsvs for each group.
#' 
#' @param mr.obj List of master regulators; list of lists, w/ positive and negative lists of MRs for each group
#' @param file.dir Directory to write file.
#' @param file.name Name to be used for each table.
#' @param num.mrs Number of MRs to write. Default of 50.
#' @param top Optional argument to write top (activated) MRs. TRUE by default.
#' @param bottom Optional argument to write bottom (deactivated) MRs. FALSE by default.
#' @export
MRTableWrite <- function(mr.obj, file.dir, file.name, num.mrs = 50, top = TRUE, bottom = FALSE) {
  file.pref <- paste(file.dir, file.name, sep = '/')
  mr.obj <- mr.obj[sort(names(mr.obj))]
  # activated mrs
  if (top) {
    active.df <- as.data.frame(lapply(mr.obj, function(x) {names(x$'positive'[1:num.mrs])} ))
    colnames(active.df) <- names(mr.obj)
    write.table(active.df, file = paste(file.pref, '_pos-mrs.tsv'), sep = '\t',
                row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  # deactivated mrs
  if (bottom) {
    deac.df <- as.data.frame(lapply(mr.obj, function(x) {names(x$'negative'[1:num.mrs])} ))
    colnames(deac.df) <- names(mr.obj)
    write.table(deac.df, file = paste(file.pref, '_pos-mrs.tsv'), sep = '\t',
                row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}


