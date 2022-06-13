#' Filters raw gene expression data based on depth, unique genes, and percentage of MT reads.
#' 
#' @param raw.counts Matrix of raw gene expression data (features X samples).
#' @param min.depth Minimum depth for each sample. Default of 0.
#' @param min.depth Minimum depth for each sample. Default of 0.
#' @param min.depth Minimum depth for each sample. Default of 0.
#' @param min.depth Minimum depth for each sample. Default of 0.
#' @param max.mt Maximum percentage of Mitochondrial reads. Default of 1.
#' @param min.gene.count Minimum total reads for a gene to be kept. If not specified, default is 0.01 * # of samples.
#' @param species One of `c('hum', 'mur')`, specifying human or murine data respectively. If not specified, assumes human.
#' @param genes One of `c('symb', 'ensg')`, specifying gene symbols or ENSG names respectively. If not specified, assumed symbols.
#' @return Matrix of normalized counts (features X samples).
#' @export
qc_filt <- function(raw.counts, min.depth = 0, max.depth = Inf, min.genes = 0, max.genes = Inf, max.mt = 1, min.gene.count,
                    species = c('hum', 'mur'), genes = c('symb', 'ensg')) {
  # check arguments
  match.arg(species)
  if (missing(species)) { species <- 'hum' }
  match.arg(genes)
  if (missing(genes)) { genes <- 'symb' }
  
  # collect statistics
  sample.depth <- colSums(raw.counts)
  sample.genes <- apply(raw.counts, 2, function(x) {length(which(x > 0))})
  mtg <- intersect(mt.genes[[paste(species, genes, sep = '.')]], rownames(raw.counts))
  sample.mt <- apply(raw.counts, 2, function(x) {
    sum(x[mtg])
  })
  sample.mt <- sample.mt / sample.depth
  qc.df <- data.frame('Depth' = sample.depth,
                      'Genes' = sample.genes,
                      'MT.per' = sample.mt,
                      'Sample' = rep('raw', length(sample.depth)))
  
  # filter cells
  keep.cells <- Reduce(intersect, list(which(sample.depth >= min.depth),
                                       which(sample.depth <= max.depth),
                                       which(sample.genes >= min.genes),
                                       which(sample.genes <= max.genes),
                                       which(sample.mt <= max.mt)))
  filt.counts <- raw.counts[, keep.cells]
  # filter genes
  if (missing(min.gene.count)) {min.gene.count <- 0.01 * ncol(filt.counts)}
  filt.counts <- filt.counts[which(rowSums(filt.counts) >= min.gene.count),]
  
  return(filt.counts)
}

#' Peforms a CPM normalization on the given matrix.
#' 
#' @param dat.mat Matrix of data, typically unnormalized counts (genes X samples).
#' @param l2 If true, will log transform the matrix. Default of FALSE.
#' @param remove.zeroes Removes rows (genes) with zero expression across all samples. Default of TRUE.
#' @return CPM normalizaed matrix.
#' @export
cpm_norm <- function(dat.mat, l2 = FALSE, remove.zeroes = TRUE) {
  # remove genes w/ zero counts
  if (remove.zeroes) { dat.mat <- dat.mat[which(rowSums(dat.mat) > 0),] }
  # log2 normalize rows
  if (l2) { dat.mat <- log2(dat.mat + 1) }
  # return matrix
  return(dat.mat)
}

#' Normalizes a matrix using one round of proportional fitting. 
#' 
#' @param dat.mat Matrix of data (features X samples).
#' @param sf Optional size factor argument. If not specified, will be computed as mean of column sums.
#' @return Normalized matrix (features X samples).
#' @export
pf_norm <- function(dat.mat, sf) {
  # comptue factors
  pf <- colSums(dat.mat)
  if (missing(sf)) {
    sf <- mean(pf)
  }
  # multiple and return matrix
  norm.mat <- apply(dat.mat, 2, function(x) {x * (sf / sum(x))})
  return(norm.mat)
}

#' Performs PF log1 PF norm normalization. 
#' 
#' @param dat.mat Matrix of data (features X samples).
#' @return Normalized matrix (features X samples).
#' @export
pflpf_norm <- function(dat.mat) {
  norm.mat <- pf_norm(dat.mat)
  norm.mat <- log(1 + norm.mat)
  norm.mat <- pf_norm(norm.mat)
  return(norm.mat)
} 

#' Generates a gene expression signature by scaling and centering, either againt a reference or internally.
#' 
#' @param dat.mat Matrix of test data (features X samples).
#' @param ref.mat Optional matrix of reference data (features X samples). If not specified, an internal GES will be generated.
#' @return GES Matrix (featues X samples).
#' @export
scale_ges <- function(dat.mat, ref.mat) {
  # get normalizing values, either from the test or reference matrix
  if (missing(ref.mat)) {
    row.means <- apply(dat.mat, 1, mean)
    row.sds <- apply(dat.mat, 1, sd)
  } else {
    # filter for shared genes
    shared.genes <- intersect(rownames(dat.mat), rownames(ref.mat))
    dat.mat <- dat.mat[shared.genes,]
    ref.mat <- ref.mat[shared.genes,]
    # get means and sd from the reference
    row.means <- apply(ref.mat, 1, mean)
    row.sds <- apply(ref.mat, 1, sd)
  }
  # transform matrix
  ges.mat <- (dat.mat - row.means) / row.sds
  return(ges.mat)
}

#' Generates a gene expression signature through an ECDF, either internally or against a reference.
#' 
#' @param dat.mat Matrix of test data (features X samples).
#' @param ref.mat Optional matrix of reference data (features X samples). If not specified, an internal GES will be generated.
#' @return GES Matrix (featues X samples).
#' @export
ecdf_ges <- function(dat.mat, ref.mat) {
  if (missing(ref.mat)) {
    ref.mat <- dat.mat
  }
  # normalize each shared gene, then reformat as matrix
  shared.genes <- intersect(rownames(dat.mat), rownames(ref.mat))
  ges.vecs <- lapply(shared.genes, function(x) {
    ecdf_norm(dat.mat[x,], ref.mat[x,])
  })
  # reformat as matrix and return
  ges.mat <- do.call(rbind, ges.vecs)
  colnames(ges.mat) <- colnames(dat.mat)
  rownames(ges.mat) <- shared.genes
  return(ges.mat)
}

#' Normalizes values in the test vector based on the ECDF of the reference vector.
#' 
#' @param test.vec Vector of test data.
#' @param ref.vec Vector of reference data.
#' @return Vector of normalized values as z-scores.
#' @export
ecdf_norm <- function(test.vec, ref.vec) {
  # fit ecdf and transform
  ecdf.func <- ecdf(ref.vec)
  ecdf.vals <- ecdf.func(test.vec)
  ecdf.step <- 1 / length(ref.vec)
  # correct 0/1 to avoid -Inf/Inf values in final GES
  ecdf.vals[which(ecdf.vals == 1)] <- 1 - ecdf.step / 2
  ecdf.vals[which(ecdf.vals == 0)] <- ecdf.step / 2
  # qnorm to get z-scores
  norm.vec <- sapply(ecdf.vals, function(y) {qnorm(y, lower.tail = TRUE)})
  return(norm.vec)
}

#' Generates a single-sample, internal GES.
#' 
#' @param filt.counts Matrix (features X samples) of filtered count data.
#' @param norm.method Normalization method. One of `c('cpm', 'pflpf')`, uses `pflpf` by default.
#' @param est.method Estimation method. One of `c('mle', 'map')`, with `map` by default.
#' @param map.iter Number of iterations to use if using MAP estimation. Default of 10.
#' @return GES Matrix (features X samples).
#' @export
internal_ges <- function(filt.counts, norm.method = c('cpm', 'pflpf'), est.method = c('map', 'ps'), map.iter = 10) {
  require(DirichletReg)
  
  # check arguments
  match.arg(norm.method)
  if (missing(norm.method)) { norm.method <- 'pflpf' }
  match.arg(est.method)
  if (missing(est.method)) { est.method <- 'map' }
  
  # add Jeffreys Prior
  filt.counts <- filt.counts + 0.5
  # mle ges
  if (est.method == 'map') {
    cat("Generating GES using an MAP...\n")
    mle.mat <- switch(norm.method, 
                      "pflpf" = pflpf_norm(filt.counts),
                      "cpm" = cpm_norm(filt.counts))
    ges.mat <- ecdf_ges(mle.mat)
  } else {
    cat("Generating a GES by sampling from the posterior. WARNING: May be slow for many samples.\n")
    map.mats <- list()
    # draw from dirichlet, then generate GES
    for (i in 1:map.iter) {
      if (i %% 10 == 0) { cat(paste("Map iteration ", i, "...\n", sep = ''))}
      dir.mat <- apply(filt.counts, 2, function(x) { rdirichlet(1, x)})
      rownames(dir.mat) <- rownames(filt.counts)
      dir.ges <- ecdf_ges(dir.mat)
      map.mats[[i]] <- dir.ges
    }
    # stouffer integrate
    ges.mat <- Reduce('+', map.mats) / sqrt(map.iter)
  }
  
  return(ges.mat)
}
