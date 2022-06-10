### PISCES SEURAT OBJECT FUNCTIONS ###
## this file contains functions designed for the Seurat object-based PISCES workflow

#' Performs louvain clustering on the active assay of the given Seurat object.
#' Looks for a distance matrix 'dist' in the 'misc' field of the active assay.
#' 
#' @param pisces.obj A Seurat object w/ either 'dist' in the 'misc' field or w/ 'scale.data'.
#' @param kmin Minimum number of neighbors Default of 5.
#' @param kmax Maximum number of neibhbors Default of 50.
#' @param kstep Step size between k values. Default of 5.
#' @return Seurat object with clustering object 'louvain.clust' added to 'misc' field of active assay.
#' And optimal clustering set as 'pisces.cluster' label for active assay.
#' @export
LouvainClust <- function(pisces.obj, kmin = 5, kmax = 50, kstep = 5) {
  # check if distance matrix has been computed
  if (!('dist' %in% names(pisces.obj@assays[[pisces.obj@active.assay]]@misc))) {
    cat("Distance matrix not found. Computing using 'cor_dist'...")
    AddDist(pisces.obj)
  }
  # cluster
  l.clust <- louvain_k(pisces.obj@assays[[pisces.obj@active.assay]]@misc$dist, 
                       kmin = kmin, kmax = kmax, kstep = kstep)
  # add to object
  pisces.obj@assays[[pisces.obj@active.assay]]@misc[['louvain.clust']] <- l.clust
  pisces.obj@meta.data[['pisces.clust']] <- l.clust$opt.clust
  return(pisces.obj)
}

#' Performs PAM clustering on the active assay of the given Seurat object.
#' Looks for a distance matrix 'dist' in the 'misc' field of the active assay.
#' 
#' @param pisces.obj A Seurat object w/ either 'dist' in the 'misc' field or w/ 'scale.data'.
#' @param kmin Minimum number of clusters. Default of 2.
#' @param kmax Maximum number of clusters. Default of 5.
#' @return Seurat object with clustering object 'louvain.clust' added to 'misc' field of active assay.
#' And optimal clustering set as 'pisces.cluster' label for active assay.
#' @export
PAMClust <- function(pisces.obj, kmin = 2, kmax = 5) {
  # check if distance matrix has been computed
  if (!('dist' %in% names(pisces.obj@assays[[pisces.obj@active.assay]]@misc))) {
    cat("Distance matrix not found. Computing using 'cor_dist'...")
    AddDist(pisces.obj)
  }
  # cluster
  p.clust <- pam_k(pisces.obj@assays[[pisces.obj@active.assay]]@misc$dist, 
                   kmin = kmin, kmax = kmax)
  # add to object
  pisces.obj@assays[[pisces.obj@active.assay]]@misc[['pam.clust']] <- p.clust
  pisces.obj@meta.data[['pisces.clust']] <- p.clust$opt.clust
  return(pisces.obj)
}

#' Adds a named assay 'PISCES' to a Seurat object.
#' 
#' @param seurat.obj Seurat object w/ RNA assay.
#' @param stage One of 'gexp' or 'pact'; adds either the 'PISCES.gexp' or 'PISCES.pact' assay. Default of 'pact'.
#' @return Seurat object w/ original assay and added 'PISCES' assay.
#' @export
AddPISCESAssay <- function(seurat.obj, stage) {
  # check for proper stage argument
  if (!(stage %in% c('gexp', 'pact'))) {
    cat("ERROR: Please set stage as either 'gexp' or 'pact'.")
    return(seurat.obj)
  }
  # add assay and set activate assay
  assay.name <- paste('PISCES', stage, sep = '')
  pisces.assay <- CreateAssayObject(counts = seurat.obj@assays$RNA@counts)
  seurat.obj[[assay.name]] <- pisces.assay
  seurat.obj@active.assay <- assay.name
  return(seurat.obj)
}

#' Adds networks to a Seurat object with a `PISCESpact` assay.
#' 
#' @param pisces.obj Seurat object with `PISCESpact`
#' @param net.list List of ARACNe3 networks.
#' @return Seurat object with `a3.nets` object added to `misc` field of `PISCESpact` assay.
#' @export
AddNetworks <- function(pisces.obj, net.list) {
  pisces.obj@assays$PISCESpact@misc[['a3.nets']] <- net.list
  return(pisces.obj)
}

#' Generates metacells from the given Seurat object. Rquires that the 'PISCESgexp' assay has been added.
#' By default, uses the 'pisces.clust' metadata to group clusters. This can be overrided with the 'clust.vec' argument.
#' If clust.vec is 'NA', then the entire matrix will be used.
#' 
#' @param pisces.obj A Seurat object with the PISCESgexp assay.
#' @param out.dir Output directory for metacell matrices to be written to.
#' @param proj.name File name preface 
#' @param clust.vec Optional clustering vector object. 
#' @param num.neighbors Number of neighbors to use in metacells. Default of 5.
#' @param subset Number of cells to subset to. Default of 250. No subsetting if set equal to NULL.
#' @param min.samps Minimum number of samples in a cluster required for meta cells. Default of 500. 
#' @param subset
#' @export
WriteMetaCells <- function(pisces.obj, file.path = '.', proj.name = '', clust.vec = NULL, num.neighbors = 5, subset = 250, min.samps = 500) {
  # check for PISCESgexp assay
  if (!('PISCESgexp' %in% names(pisces.obj@assays))) {
    cat("Error: No 'PISCESgexp' assay detected. Please see 'AddPISCESAssay'.")
  }
  # establish clustering vector
  if (is.null(clust.vec)) {
    cat("Using 'pisces.clust' as clustering vector.\n")
    clust.vec <- pisces.obj$pisces.clust
  } else if (is.na(clust.vec)) {
    cat("Generating one metacell matrix from all cells.\n")
    clust.vec <- rep('all', ncol(pisces.obj@assays$PISCESgexp@counts))
    names(clust.vec) <- colnames(pisces.obj@assays$PISCESgexp@counts)
  } 
  # generate metacells
  meta.mats <- make_metacells(pisces.obj@assays$PISCESgexp@counts, pisces.obj@assays$PISCESgexp@misc$dist,
                              clust.vec, num.neighbors, subset, min.samps)
  # write to specified file path
  cat("Writing to file...")
  for (m.name in names(meta.mats)) {
    saveRDS(meta.mats[[m.name]], file = paste(file.path, '/', proj.name, '_', m.name, '-meta-cells.rds', sep = ''))
  }
}

#' Generates protein activity readout using the NaRnEA algorithm.
#' Uses either the GES object found in the `misc` field of the `PISCESpact` assay or
#' the `scale.data` field of the `PISCESgexp` assay.
#' 
#' @param pisces.obj Seurat object w/ the `PISCESpact` assay. 
#' @param net.list List of networks. If not specified, networks in the `misc` field of the `PISCESpact` assay will be used.
#' @param sample.weights Flag to compute sample-specific weights. Default of TRUE.
#' @param use.sct Flag to use the SCT matrix as a GES. Default of FALSE.
#' @return Seurat object with `PES` matrix in `scale.data` of `PISCESpact` assay.
#' `NES` results will be stored in the `misc` field of the `PISCESpact` assay.
#' @export
PISCESActivity <- function(pisces.obj, net.list, sample.weights = TRUE, use.sct = FALSE) {
  # identify GES matrix
  if (use.sct) {
    ges.mat <- pisces.obj@assays$SCT@scale.data
  } else if ('ges' %in% names(pisces.obj@assays$PISCESpact@misc)) {
    ges.mat <- pisces.obj@assays$PISCESpact@misc$ges
  } else {
    ges.mat <- pisces.obj@assays$PISCESgexp@scale.data
  }
  # identify networks
  if (missing(net.list)) {
    if ('a3.nets' %in% names(pisces.obj@assays$PISCESpact@misc)) {
      net.list <- pisces.obj@assays$PISCESpact@misc
    } else {
      cat("No networks found in `pisces.obj` - please add or specify in the call to `PISCESActivtiy`\n")
      return(pisces.obj)
    }
  }
  # run narnea
  narnea.res <- meta_narnea(ges.mat, net.list, sample.weights)
  # add to object
  pisces.obj@assays$PISCESpact@scale.data <- narnea.res$PES
  pisces.obj@assays$PISCESpact@misc[['NES']] <- narnea.res$NES
  return(pisces.obj)
}

#' Generates a principle component analysis (PCA) for the active assay of the given object.
#' 
#' @param pisces.obj A Seurat object with 'scale.data' in the active assay.
#' @param num.pcs Flag to use principle components. Default of 100.
#' @return A Seurat object with a 'pca' object added to the 'misc' field of the active assay.
#' @export
AddPCA <- function(pisces.obj, num.pcs = 100) {
  # generate pca
  pca.obj <- fast_pca(pisces.obj@assays[[pisces.obj@active.assay]]@scale.data, num.pcs)
  # add to object and return
  pisces.obj@assays[[pisces.obj@active.assay]]@misc[['pca']] <- pca.obj
  return(pisces.obj)
}

#' Generates a Multidimensional Scaling (MDS) for the active assay of the given object.
#' Checks for a distance object in the 'misc' field; calculates one based on correlation distance if not present.
#' 
#' @param pisces.obj A Seurat object with 'scale.data' in the active assay.
#' @param num.dims Number of dimensions to return. Default of 2.
#' @return A Seurat object with a 'mds' object added to the 'misc' field of the active assay.
#' @export
AddMDS <- function(pisces.obj, num.dims = 2) {
  # check for distance object
  if (!('dist' %in% names(pisces.obj@assays[[pisces.obj@active.assay]]@misc))) {
    cat("No distance matrix found - calculating new distance matrix...\n")
    dist.mat <- cor_dist(pisces.obj@assays[[pisces.obj@active.assay]]@scale.data)
  } else {
    dist.mat <- pisces.obj@assays[[pisces.obj@active.assay]]@misc$dist
  }
  # generate mds
  mds.mat <- cmdscale(dist.mat, k = num.dims)
  # add to object
  pisces.obj@assays[[pisces.obj@active.assay]]@misc[['mds']] <- mds.mat
  return(pisces.obj)
}

#' Generates a UMAP for the active assay of the given object.
#' 
#' @param pisces.obj A Seurat object with 'scale.data' in the active assay.
#' @param dense.map Flag to return dense-UMAP. Default of FALSE.
#' @return A Seurat object with a 'umap' object added to the 'misc' field of the active assay.
#' @export
AddUMAP <- function(pisces.obj, dense.map = FALSE) {
  # generate umap
  if (dense.map) {
    umap.mat <- make_umap(pisces.obj@assays[[pisces.obj@active.assay]]@scale.data)
  } else {
    umap.mat <- make_umap(pisces.obj@assays[[pisces.obj@active.assay]]@scale.data)
  }
  # add to object
  pisces.obj@assays[[pisces.obj@active.assay]]@misc[['umap']] <- umap.mat
  return(pisces.obj)
}

#' Generates a distance matrix for the active assay of the given object.
#' 
#' @param pisces.obj A Seurat object w/ either 'scale.data' in the active assay or 'pca' object in the 'misc' field of the active assay.
#' @param use.pcs Flag to indicate whether PCA should be used. If true, object must have 'pca' object in the 'misc' field of the active assay.
#' @param dist.method One of 'cor', 'l1', or 'l2'. Default of 'cor'.
#' @return Seurat object with a 'dist' object added to the 'misc' field of the active assay.
#' @export
AddDist <- function(pisces.obj, use.pcs = TRUE, dist.method = c('cor', 'l1', 'l2')) {
  match.arg(dist.method)
  if (missing(dist.method)) { dist.method <- 'cor' }
  
  # identify data matrix
  if (use.pcs) {
    if (!('pca' %in% names(pisces.obj@assays[[pisces.obj@active.assay]]@misc))) {
      cat("Error: use.pcs is TRUE, but no 'pca' field in 'misc' of active assay.\n")
      return(pisces.obj)
    }
    dat.mat <- pisces.obj@assays[[pisces.obj@active.assay]]@misc$pca$x
  } else {
    dat.mat <- pisces.obj@assays[[pisces.obj@active.assay]]@scale.data
  }
  
  # generate distance matrix
  dist.mat <- switch (dist.method,
    'cor' = cor_dist(t(dat.mat)),
    'l1' = dist(dat.mat, method = 'manhattan'),
    'l2' = dist(dat.mat, method = 'euclidean')
  )
  
  # add to object and return
  pisces.obj@assays[[pisces.obj@active.assay]]@misc[['dist']] <- dist.mat
  return(pisces.obj)
}

#' Normalizes gene expression in the PISCESgexp assay. Stores normalized counts in 'data' field.
#' Additionally, generates a scaled matrix (interpretable as an internal GES) to the 'scale.data' field.
#' 
#' @param pisces.obj A Seurat object with the PISCESgexp assay.
#' @param method One of 'cpm' or 'pflogpf'
#' @param filt.per Threshold value for gene filtration; genes w/ less than filt.per * num.samps total counts will be removed. Default of 0.01.
#' @return A Seurat object with normalized expression in the 'data' field of the PISCESgexp assay.
#' @export
NormalizeExpression <- function(pisces.obj, method = c('cpm', 'pflogpf'), filt.per = 0.01) {
  # check method argument
  match.arg(method)
  if (missing(method)) { method <- 'pflogpf' }
  # check that pisces.obj has the PISCESgexp assay
  if (!('PISCESgexp' %in% names(pisces.obj@assays))) {
    cat("Error: seurat object must contain the 'PISCESgexp' assay\n")
    return(pisces.obj)
  }
  # perform normalization
  counts.mat <- as.matrix(pisces.obj@assays$PISCESgexp@counts)
  count.thresh <- floor(ncol(counts.mat) * filt.per) - 1
  cat(paste("Removing genes with fewer than", count.thresh, "total reads...\n"))
  counts.mat <- counts.mat[which(rowSums(counts.mat) > count.thresh),]
  pisces.obj@assays$PISCESgexp@counts <- counts.mat
  if (method == 'cpm') {
    pisces.obj@assays$PISCESgexp@data <- cpm_norm(counts.mat, l2 = TRUE)
  } else if (method == 'pflogpf') {
    pisces.obj@assays$PISCESgexp@data <- pflpf_norm(counts.mat)
  }
  pisces.obj@assays$PISCESgexp@scale.data <- scale_ges(pisces.obj@assays$PISCESgexp@data)
  return(pisces.obj)
}

#' Generates a GES, either internally or against an external reference.
#' 
#' @param pisces.obj A Seurat object with the `PISCESpact` assay.
#' @param method One of `scale` or `ecdf`. If not specified, will use `scale`.
#' @param ref.mat Optional matrix of appropriately normalized reference data (features X samples).
#' @return A Seurat object with a `ges` matrix (features X samples) added to the `PISCESpact` assay.
#' @export
AddGES <- function(pisces.obj, method = c('scale', 'ecdf'), ref.mat = NULL) {
  # check method argument
  match.arg(method)
  if (missing(method)) { method <- 'scale' }
  # add appropriate GES
  if (is.null(ref.mat)) {
    if (method == 'scale') { ges.mat <- pisces.obj@assays$PISCESgexp@scale.data }
    if (method == 'ecdf') { ges.mat <- ecdf_ges(pisces.obj@assays$PISCESgexp@counts) }
  } else {
    if (method == 'scale') { ges.mat <- scale_ges(pisces.obj@assays$PISCESgexp@counts,
                                                  ref.mat)}
    if (method == 'ecdf') { ges.mat <- ecdf_ges(pisces.obj@assays$PISCESgexp@counts,
                                                ref.mat)}
  }
  pisces.obj@assays$PISCESpact@misc[['ges']] <- ges.mat
  return(pisces.obj)
}
