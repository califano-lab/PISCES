# Shared fixtures for the PISCES test suite.
#
# testthat sources every `helper-*.r` file before running the tests, so the
# functions below are available in all test files.

suppressPackageStartupMessages({
  library(Seurat)
  library(R6)
})

# The classes source each other with paths relative to the repository root,
# while testthat runs tests from tests/testthat, so load them from the root.
pisces_root <- normalizePath(file.path("..", ".."))
local({
  owd <- setwd(pisces_root)
  on.exit(setwd(owd))
  suppressPackageStartupMessages({
    source("classes/metacell_generator.r")
    source("classes/runner.r")
  })
})

#' Simulate a Small Clustered Seurat Object with Known Regulatory Structure
#'
#' Each cluster gets its own baseline expression. Within a cluster, every
#' regulator has a latent per-cell activity, and its targets follow that
#' activity, so ARACNe3 has genuine regulator -> target dependencies to find.
#' Counts are Poisson, and per-cell sequencing depth varies, as in real data.
#'
#' The object carries what the pipeline reads:
#' - an "SCT" assay with a counts layer (MetacellGenerator),
#' - `seurat_clusters` metadata (MetacellGenerator),
#' - an "integrated" assay with a scale.data layer (Runner$run_viper()).
#'
#' @param cluster_sizes   Number of cells per cluster.
#' @param n_regulators    Number of regulators.
#' @param targets_per_reg Number of targets driven by each regulator.
#' @param n_background    Number of genes with no regulator.
#' @param seed            Seed for the simulation.
#'
#' @return A list with `seurat_obj` and `regulators` (character vector).
#' Regulators get 30 targets by default because VIPER ignores regulons with
#' fewer than 25 targets (its `minsize` default).
simulate_seurat <- function(cluster_sizes = c(500, 300, 80),
                            n_regulators = 12, targets_per_reg = 30,
                            n_background = 40, seed = 1) {
  set.seed(seed)
  regulators <- sprintf("REG%02d", seq_len(n_regulators))
  targets <- sprintf("TGT%03d", seq_len(n_regulators * targets_per_reg))
  background <- sprintf("BKG%03d", seq_len(n_background))
  genes <- c(regulators, targets, background)
  target_owner <- rep(seq_len(n_regulators), each = targets_per_reg)

  blocks <- lapply(seq_along(cluster_sizes), function(k) {
    n <- cluster_sizes[k]
    base <- rnorm(length(genes), mean = 1.5, sd = 0.7)
    activity <- matrix(rnorm(n_regulators * n), nrow = n_regulators)
    log_mu <- matrix(base, nrow = length(genes), ncol = n)
    log_mu[seq_len(n_regulators), ] <-
      log_mu[seq_len(n_regulators), ] + 0.8 * activity
    target_rows <- n_regulators + seq_along(targets)
    log_mu[target_rows, ] <-
      log_mu[target_rows, ] + 0.8 * activity[target_owner, , drop = FALSE]
    depth <- rlnorm(n, meanlog = 0, sdlog = 0.3)
    mu <- exp(log_mu) * rep(depth, each = length(genes))
    matrix(rpois(length(mu), mu), nrow = length(genes))
  })

  counts <- do.call(cbind, blocks)
  dimnames(counts) <- list(genes, sprintf("cell%04d", seq_len(ncol(counts))))
  clusters <- factor(rep(seq_along(cluster_sizes) - 1, cluster_sizes))

  obj <- suppressWarnings(CreateSeuratObject(counts = counts))
  obj[["SCT"]] <- CreateAssayObject(counts = counts)

  log_norm <- log1p(t(t(counts) / colSums(counts)) * 1e4)
  obj[["integrated"]] <- CreateAssayObject(data = log_norm)
  obj <- SetAssayData(obj, assay = "integrated", layer = "scale.data",
                      new.data = t(scale(t(log_norm))))
  obj$seurat_clusters <- unname(clusters)

  list(seurat_obj = obj, regulators = regulators,
       planted_edges = data.frame(regulator = regulators[target_owner],
                                  target = targets))
}

#' Run the Metacell Step into a Fresh Temporary Directory
#'
#' @param sim  Output of `simulate_seurat()`.
#' @param ...  Passed to `generate_metacell_matrices()`.
#'
#' @return A list with the returned matrices (`mats`) and the directory
#'         (`dir`).
make_metacells <- function(sim, ...) {
  dir <- tempfile("metacells_")
  dir.create(dir)
  gen <- MetacellGenerator$new(sim$seurat_obj, dir)
  list(mats = gen$generate_metacell_matrices(...), dir = dir)
}

#' Read a Metacell TSV Written by MetacellGenerator Back into a Matrix
read_metacell_tsv <- function(path) {
  as.matrix(read.table(path, header = TRUE, sep = "\t", row.names = 1,
                       check.names = FALSE))
}

#' Path to a Real ARACNe3 Binary, or "" if None Is Configured
#'
#' Integration tests run only when the environment variable ARACNE3_BIN points
#' to an executable ARACNe3 build.
aracne3_binary <- function() {
  bin <- Sys.getenv("ARACNE3_BIN", unset = "")
  if (nzchar(bin) && file.exists(bin)) bin else ""
}
