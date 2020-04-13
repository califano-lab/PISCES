#' Transforms gene names from Ensembl to hgnc.
#' 
#' @param dat.mat Matrix of data with ENSEMBL names (genes X samples).
#' @return Data with HGNC names. Some data will likely be lost in the conversion.
#' @export
Ensemble2GeneName<-function(dat.mat) {
  # get the mart
  ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "http://useast.ensembl.org")
  name.map <- biomaRt::getBM(attributes=c('hgnc_symbol','hgnc_id','ensembl_gene_id'), 
                    filters = 'ensembl_gene_id', values = (rownames(dat.mat)), mart = ensembl)
  # remove rows with no match for gene name, then match and replace
  name.map <- name.map[ which(!is.na(name.map$hgnc_symbol)) , ]
  name.map <- name.map[ which(name.map$hgnc_symbol != '') , ]
  convert.dat <- dat.mat[ name.map$ensembl_gene_id , ]
  rownames(convert.dat) <- name.map$hgnc_symbol
  return(convert.dat)
}

#' Transforms gene names from Ensembl to Entrez
#' 
#' @param dat.mat Matrix of data with ENSEMBL names (genes X samples).
#' @return Data with Entrez names. Some data will likely be lost in the conversion.
#' @export
Ensemble2Entrez <- function(dat.mat) {
  # get the mart
  ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "http://useast.ensembl.org")
  name.map <- biomaRt::getBM(attributes=c('entrezgene_id','ensembl_gene_id'), 
                    filters = 'ensembl_gene_id', values = (rownames(dat.mat)), mart = ensembl)
  # remove rows with no match for entrez, then match and replace
  name.map <- name.map[ which(!is.na(name.map$entrezgene)) , ]
  name.map <- name.map[ which(name.map$entrezgene != '') , ]
  convert.dat <- dat.mat[ name.map$ensembl_gene_id , ]
  rownames(convert.dat) <- name.map$entrezgene
  return(convert.dat)
}

#' Transforms gene names from Entrez to Ensembl
#' 
#' @param dat.mat Matrix of data with Entrez names (genes X samples).
#' @return Data with Ensembl names. Some data will likely be lost in the conversion.
#' @export
Entrez2Ensemble <- function(dat.mat) {
  # get the mart
  ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "http://useast.ensembl.org")
  name.map <- biomaRt::getBM(attributes=c('entrezgene_id','ensembl_gene_id'), 
                    filters = 'entrezgene_id', values = (rownames(dat.mat)), mart = ensembl)
  # remove rows with no match for entrez, then match and replace
  name.map <- name.map[ which(!is.na(name.map$ensembl_gene_id)) , ]
  name.map <- name.map[ which(name.map$ensembl_gene_id != '') , ]
  convert.dat <- dat.mat[ as.character( name.map$entrezgene ) , ]
  rownames(convert.dat) <- name.map$ensembl_gene_id
  return(convert.dat)
}