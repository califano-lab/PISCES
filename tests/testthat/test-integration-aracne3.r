# End-to-end tests with a real ARACNe3 build. They run only when the
# environment variable ARACNE3_BIN points to an ARACNe3 executable, e.g.
#   ARACNE3_BIN=/path/to/ARACNe3_app_release Rscript tests/run_tests.r

sim <- simulate_seurat()
metacells <- make_metacells(sim)

run_real_pipeline <- function(num_subnets = 1) {
  root <- tempfile("pipeline_")
  dirs <- file.path(root, c("out", "out/aracne_results", "out/viper_results"))
  for (d in dirs) dir.create(d, recursive = TRUE)
  reg_file <- file.path(root, "tfs.txt")
  writeLines(sim$regulators, reg_file)
  runner <- Runner$new(metacells$mats, dirs[1], aracne3_binary(), dirs[2],
                       list(tfs = reg_file), sim$seurat_obj, dirs[3],
                       threads = 2, seed = 42, num_subnets = num_subnets)
  # ARACNe3 caches its null model in the working directory.
  withr::with_dir(root, capture.output(runner$run_aracne()))
  runner
}

read_network <- function(runner, mat_name) {
  read.table(file.path(runner$aracne_output_path, paste0("tfs_", mat_name),
                       "consolidated-net_defaultid.tsv"),
             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

test_that("ARACNe3 receives the subsampled metacells, not all of them", {
  skip_if(aracne3_binary() == "", "ARACNE3_BIN is not set")
  runner <- run_real_pipeline()

  for (mat_name in names(metacells$mats)) {
    log <- readLines(file.path(runner$aracne_output_path,
                               paste0("tfs_", mat_name), "subnets_log",
                               "log_subnet1_defaultid.txt"))
    n_samples <- as.integer(sub(".*: ", "",
                                grep("^Total # samples", log, value = TRUE)))
    expect_equal(n_samples, ncol(metacells$mats[[mat_name]]))
  }
})

test_that("ARACNe3 recovers the planted regulator -> target edges", {
  skip_if(aracne3_binary() == "", "ARACNE3_BIN is not set")
  runner <- run_real_pipeline()
  planted <- paste(sim$planted_edges$regulator, sim$planted_edges$target)

  # Precision is deliberately not held near 1: kNN metacells overlap, which
  # adds false edges (in this fixture about 80% precision for the 500-cell
  # cluster and 60% for the 300-cell one, where 200 of 300 metacells are
  # kept; 200 single cells give about 93%). A random edge is planted with
  # probability of about 7%, so 50% still shows the networks are driven by
  # the planted structure.
  for (mat_name in names(metacells$mats)[1:2]) {  # the 200-metacell clusters
    net <- read_network(runner, mat_name)
    found <- paste(net$regulator.values, net$target.values)
    expect_gt(mean(planted %in% found), 0.9)   # recall
    expect_gt(mean(found %in% planted), 0.5)   # precision
  }
})

test_that("num_subnets sets how many subnetworks ARACNe3 builds", {
  skip_if(aracne3_binary() == "", "ARACNE3_BIN is not set")
  runner <- run_real_pipeline(num_subnets = 3)

  run_dir <- file.path(runner$aracne_output_path,
                       paste0("tfs_", names(metacells$mats)[1]))
  expect_length(list.files(file.path(run_dir, "subnets")), 3)
  # With several subnetworks, edges get different p-values.
  net <- read_network(runner, names(metacells$mats)[1])
  expect_gt(length(unique(net$log.p.values)), 1)
})

test_that("the full ARACNe -> VIPER pipeline runs", {
  skip_if(aracne3_binary() == "", "ARACNE3_BIN is not set")
  skip_if_not_installed("viper")
  runner <- run_real_pipeline()
  suppressMessages(suppressWarnings(capture.output(runner$run_viper())))

  viper_res <- readRDS(file.path(runner$viper_output_path,
                                 "viper_results.rds"))
  expect_equal(colnames(viper_res), colnames(sim$seurat_obj))
  expect_true(all(rownames(viper_res) %in% sim$regulators))
  expect_true(all(is.finite(viper_res)))
})
