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

#' Converts a list of gene names from 'from' to 'to' conventions.
#' 
#' @param gene.list List of genes to be converted.
#' @param species Species for this data. Either 'human' or 'mouse'.
#' @param from Starting name convention. One of 'ensg', 'gn', or 'entrez'.
#' @param to Ending name convention. One of 'ensg', 'gn', or 'entrez'.
#' @return Matrix with converted row names.
#' @export
convert_genes <- function(gene.list, species, from, to) {
  require(biomaRt)
  
  attr.list <- list('ensg' = 'ensembl_gene_id', 'gn' = 'hgnc_symbol', 'entrez' = 'entrezgene_id')
  # get mart for appropriate species
  if (species == 'human') {
    ensembl.mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  } else if (species == 'mouse') {
    ensembl.mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    attr.list['gn'] <- 'mgi_symbol'
  } else {
    print(paste('Species', species, 'not supported. Exiting'))
    return()
  }
  # get name.map based on from and to
  attr.set <- unlist(attr.list[c(from, to)])
  name.map <- biomaRt::getBM(attributes = attr.set, filters = attr.set[1],
                             values = gene.list, mart = ensembl.mart)
  # filter and return
  ret.list <- name.map[,2]
  ret.list <- ret.list[ which(ret.list != '') ]
  ret.list <- ret.list[ which(!is.na(ret.list)) ]
  return(ret.list)
}

#' Converts the gene rownames of a matrix from the `from` to `to` convention.
#' Removes rows of the matrix for which no name in the `to` convention can be found.
#'
#' @param dat.mat Matrix with rownames to be changed. 
#' @param species Species for this data. Either 'human' or 'mouse'.
#' @param from Starting name convention. One of 'ensg', 'gn', or 'entrez'.
#' @param to Ending name convention. One of 'ensg', 'gn', or 'entrez'.
#' @return Matrix with converted row names.
#' @export
matrix_convert_genes <- function(dat.mat, species, from, to) {
  attr.list <- list('ensg' = 'ensembl_gene_id', 'gn' = 'hgnc_symbol', 'entrez' = 'entrezgene_id')
  # get mart for appropriate species
  if (species == 'human') {
    ensembl.mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  } else if (species == 'mouse') {
    ensembl.mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    attr.list['gn'] <- 'mgi_symbol'
  } else {
    print(paste('Species', species, 'not supported. Exiting'))
    return()
  }
  # get name.map based on from and to
  attr.set <- unlist(attr.list[c(from, to)])
  name.map <- biomaRt::getBM(attributes = attr.set, filters = attr.set[1],
                             values = rownames(dat.mat), mart = ensembl.mart)
  # replace names and return
  name.map <- name.map[ which(!is.na(name.map[,2])) ,]
  name.map <- name.map[ which(name.map[,2] != '') ,]
  name.map <- name.map[ match(unique(name.map[,2]), name.map[,2]) ,]
  convert.dat <- dat.mat[ as.character(name.map[,1]) ,]; rownames(convert.dat) <- name.map[,2]
  return(convert.dat)
}
