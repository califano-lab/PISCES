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
    s.int <- sum(dat.vect) / sqrt(length(dat.vect))
  }
  return(s.int)
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
