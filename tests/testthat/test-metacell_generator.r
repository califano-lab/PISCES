sim <- simulate_seurat()  # clusters of 500, 300 and 80 cells
cluster_sizes <- c(500, 300, 80)

test_that("returns one named, subsampled matrix per cluster", {
  res <- make_metacells(sim)

  expect_named(res$mats, paste0("metacell_clust-", 1:3))
  # 200 metacells for the large clusters; a cluster with fewer cells than
  # sub_size keeps all of its metacells.
  expect_equal(unname(vapply(res$mats, ncol, integer(1))), c(200L, 200L, 80L))
})

test_that("returned matrices are CPM-normalised", {
  res <- make_metacells(sim)
  for (mat in res$mats) {
    expect_equal(unname(colSums(mat)), rep(1e6, ncol(mat)))
  }
})

test_that("writes exactly the _all and _sub files, with clean names", {
  res <- make_metacells(sim)
  expected <- unlist(lapply(1:3, function(i) {
    paste0("metacell_clust-", i, "-metaCells_",
           rep(c("all", "sub"), each = 2), c(".rds", ".tsv"))
  }))
  expect_setequal(list.files(res$dir), expected)
})

test_that("the _sub files hold exactly the returned matrix", {
  res <- make_metacells(sim)
  for (i in 1:3) {
    prefix <- file.path(res$dir, paste0("metacell_clust-", i, "-metaCells_sub"))
    returned <- res$mats[[i]]
    expect_identical(readRDS(paste0(prefix, ".rds")), returned)
    expect_equal(read_metacell_tsv(paste0(prefix, ".tsv")), returned,
                 tolerance = 1e-12)
  }
})

test_that("the _all files hold every metacell as raw counts", {
  res <- make_metacells(sim)
  for (i in 1:3) {
    all_mat <- readRDS(file.path(res$dir, paste0("metacell_clust-", i,
                                                 "-metaCells_all.rds")))
    expect_equal(ncol(all_mat), cluster_sizes[i])
    expect_true(all(all_mat == round(all_mat)))  # summed counts, not CPM

    # Each returned metacell is the CPM of the same metacell in _all.
    sub_mat <- res$mats[[i]]
    expect_true(all(colnames(sub_mat) %in% colnames(all_mat)))
    from_all <- all_mat[, colnames(sub_mat)]
    expect_equal(sub_mat, t(t(from_all) / colSums(from_all) * 1e6))
  }
})

test_that("each metacell sums the cell and its num_neighbors neighbours", {
  res <- make_metacells(sim, num_neighbors = 5)
  all_mat <- readRDS(file.path(res$dir, "metacell_clust-1-metaCells_all.rds"))

  # Recompute independently: neighbours by 1 - Pearson correlation within
  # the cluster, then the cell plus its 5 nearest neighbours summed.
  counts <- as.matrix(sim$seurat_obj[["SCT"]]@counts)[, 1:500]
  counts <- counts[rownames(all_mat), ]
  dist_mat <- 1 - cor(counts)
  expected <- vapply(seq_len(ncol(counts)), function(j) {
    neighbours <- order(dist_mat[j, ])[2:6]
    rowSums(counts[, c(j, neighbours)])
  }, numeric(nrow(counts)))
  dimnames(expected) <- dimnames(counts)

  expect_equal(all_mat, expected)
})

test_that("the subsample is reproducible and controlled by seed", {
  first <- make_metacells(sim)$mats
  again <- make_metacells(sim)$mats
  other <- make_metacells(sim, seed = 7)$mats

  expect_identical(first, again)
  expect_false(identical(colnames(first[[1]]), colnames(other[[1]])))
  # The 80-cell cluster is not subsampled, so it does not depend on the seed.
  expect_setequal(colnames(first[[3]]), colnames(other[[3]]))
})

test_that("seeding leaves the caller's random number stream untouched", {
  set.seed(123)
  expected <- runif(3)

  set.seed(123)
  make_metacells(sim)
  expect_identical(runif(3), expected)
})

test_that("seeding does not create a seed if the session had none", {
  global <- globalenv()
  saved <- get(".Random.seed", envir = global)
  on.exit(assign(".Random.seed", saved, envir = global))

  rm(".Random.seed", envir = global)
  make_metacells(sim)
  expect_false(exists(".Random.seed", envir = global, inherits = FALSE))
})

test_that("seed = NULL draws from the session's random number stream", {
  set.seed(1)
  a <- make_metacells(sim, seed = NULL)$mats
  b <- make_metacells(sim, seed = NULL)$mats
  set.seed(1)
  a_again <- make_metacells(sim, seed = NULL)$mats

  expect_false(identical(colnames(a[[1]]), colnames(b[[1]])))
  expect_identical(a, a_again)
})

test_that("sub_size above 500 is honoured and .rds and .tsv agree", {
  # Regression test: the files used to be cut to 500 columns, sampled
  # separately for the .rds and the .tsv.
  big <- simulate_seurat(cluster_sizes = c(650, 60), seed = 2)
  res <- make_metacells(big, sub_size = 600)
  prefix <- file.path(res$dir, "metacell_clust-1-metaCells_sub")

  expect_equal(ncol(res$mats[[1]]), 600)
  expect_identical(readRDS(paste0(prefix, ".rds")), res$mats[[1]])
  expect_equal(colnames(read_metacell_tsv(paste0(prefix, ".tsv"))),
               colnames(res$mats[[1]]))
})

test_that("clusters at or below size_thresh are skipped", {
  res <- make_metacells(sim, size_thresh = 80)
  expect_named(res$mats, paste0("metacell_clust-", 1:2))
})

test_that("invalid sub_size and seed are rejected", {
  gen <- MetacellGenerator$new(sim$seurat_obj, tempfile())
  expect_error(gen$generate_metacell_matrices(sub_size = 0), "sub_size")
  expect_error(gen$generate_metacell_matrices(sub_size = 10.5), "sub_size")
  expect_error(gen$generate_metacell_matrices(sub_size = c(1, 2)), "sub_size")
  expect_error(gen$generate_metacell_matrices(seed = -1), "seed")
})
