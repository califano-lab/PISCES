#' Runs VIPER on a seurat.obj that has the appropriately generated 'PISCES' assay and scaled data.
#' 
#' @param data.object Seurat object w/ 'PISCES' assay and scaled data.
#' @param net.list List of networks OR a single network.
#' @param sct.ges Optional flag to use SCT scale data as GES. Default of False.
#' @return Seurat.object with added 'viper' matrix in PISCES assay
#' @export
PISCESViper <- function(data.object, net.list, sct.ges = FALSE) {
  # check for PISCES assay
  if (!HasPISCESAssay(data.object)) {
    print('Error: No PISCES assay detected')
    return(data.object)
  }
  # get GES
  if (sct.ges) {
    ges.mat <- as.matrix(data.object@assays$SCT@scale.data)
  } else {
    ges.mat <- data.object@assays$PISCES@misc$GES
  }
  # run viper
  vip.mat <- viper::viper(ges.mat, net.list, 
                          method = 'none', eset.filter = FALSE)
  # add to object
  data.object@assays$PISCES@scale.data <- vip.mat
  return(data.object)
}

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

#' MetaVIPER implementation that will perform a weighted stouffer integration based on highest NES.
#' 
#' @param ges Gene Expression Signature (features X samples)
#' @param net.list List object with the networks to be used
#' @param use.nets Optional argument to sslect the top n networks. If not specified, all networks are used.  
#' @param ret.weights Optional argument to return the network weight matrix as well as the VIPER matrix. FALSE by default.
#' @return Either a viper matrix, or a list with a viper matrix and the network weight matrix.
WeightedVIPER <- function(ges, net.list, use.nets, ret.weights = FALSE) {
  require(viper)
  num.nets <- length(net.list)
  num.samps <- ncol(ges)
  ## create weight matrix
  w.mat <- matrix(0L, nrow = num.nets, ncol = ncol(ges))
  colnames(w.mat) <- colnames(ges); rownames(w.mat) <- names(net.list)
  ## run VIPER with each network
  print('Generating VIPER matrices...')
  vip.list <- list()
  for (i in 1:num.nets) {
    vip.list[[i]] <- viper(ges, net.list[i], method = 'none')
  }
  names(vip.list) <- names(net.list)
  ## count for each gene
  print('Generating weights...')
  uni.genes <- unique(unlist(lapply(vip.list, rownames)))
  for (g in uni.genes) {
    for (s in 1:num.samps) {
      nes.vals <- unlist(lapply(vip.list, function(x){ 
        if (g %in% rownames(x)) {
          return(x[g,s])
        } else {
          return(0)
        }}))
      max.ind <- which.max(abs(nes.vals))
      w.mat[max.ind, s] <- w.mat[max.ind, s] + 1
    }
  }
  ## integration
  print('Integrating...')
  int.mat <- matrix(0L, nrow = length(uni.genes), ncol = num.samps)
  rownames(int.mat) <- uni.genes; colnames(int.mat) <- colnames(ges)
  for (g in uni.genes) {
    for (s in 1:num.samps) {
      nes.vals <- unlist(lapply(vip.list, function(x){ 
        if (g %in% rownames(x)) {
          return(x[g,s])
        } else {
          return(NA)
        }}))
      w.vals <- w.mat[,s][!is.na(nes.vals)]
      w.vals <- w.vals / sum(w.vals)
      nes.vals <- nes.vals[!is.na(nes.vals)]
      # if use.nets are specified, subset to the top n (or use all if n > length)
      if (!missing(use.nets)) {
        w.order <- order(w.vals, decreasing = TRUE)
        w.vals <- w.vals[ w.order[1:min(length(w.order), use.nets)] ]
        nes.vals <- nes.vals[ w.order[1:min(length(nes.vals), use.nets)] ]
      }
      int.mat[g,s] <- sum(nes.vals * w.vals) / sqrt(sum(w.vals**2))
    }
  }
  ## return
  if (ret.weights) {
    return( list('viper' = int.mat, 'weights' = w.mat) )
  } else {
    return( int.mat )
  }
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