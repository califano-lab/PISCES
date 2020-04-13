#' Read in data in 10x format
#'
#' @param dat.path Path to directory containing 'matrix.mtx', 'genes.tsv', and 'barcodes.tsv'.
#' @return Matrix of raw gene expression (genes X samples).
#' @export
read10X <- function(dat.path) {
  # read in data
  raw.mat <- Matrix::as.matrix(Matrix::readMM( paste(dat.path, 'matrix.mtx', sep = '') ))
  genes <- read.table( paste(dat.path, 'genes.tsv', sep = ''), sep = '\t')
  barcodes <- read.table( paste(dat.path, 'barcodes.tsv', sep = ''), sep = '\t')
  # set names and return
  colnames(raw.mat) <- barcodes[,1]
  rownames(raw.mat) <- genes[,1]
  return(raw.mat)
}

#' Saves a matrix in a format for input to ARACNe
#'
#' @param dat.mat Matrix of data (genes X samples).
#' @param out.file Output file where matrix will be saved.
#' @param subset Switch for subsetting the matrix to 500 samples. Default TRUE.
#' @export
ARACNeTable <- function(dat.mat, out.file, subset = TRUE) {
  dat.mat <- dat.mat[!duplicated(rownames(dat.mat)), ]
  if (subset) {
    dat.mat <- dat.mat[, sample(colnames(dat.mat), min(ncol(dat.mat), 500)) ]
  }
  sample.names <- colnames(dat.mat)
  gene.ids <- rownames(dat.mat)
  m <- dat.mat
  mm <- rbind( c("gene", sample.names), cbind(gene.ids, m))
  write.table( x = mm , file = out.file ,
               sep="\t", quote = F , row.names = F , col.names = F )
}
