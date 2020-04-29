#' Converts gene names from 'from' to 'to' convenctions.
#'
#' @param dat.mat Matrix with rownames to be changed.
#' @param species Species for this data. Either 'human' or 'mouse'.
#' @param from Starting name convention. One of 'ensg', 'gn', or 'entrez'.
#' @param to Ending name convention. One of 'ensg', 'gn', or 'entrez'.
#' @return Matrix with converted row names.
#' @export
GeneNameConvert <- function(dat.mat, species, from, to) {
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
  convert.dat <- dat.mat[ name.map[,1] ,]; rownames(convert.dat) <- name.map[,2]
  return(convert.dat)
}
