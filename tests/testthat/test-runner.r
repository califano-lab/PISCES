sim <- simulate_seurat()
metacells <- make_metacells(sim)

#' Create a Stand-in ARACNe3 Executable
#'
#' The script appends its arguments to `calls.log` next to it, and writes a
#' consolidated network (the planted regulator -> target edges) into the `-o`
#' directory, as ARACNe3 does. With `exit_status` != 0 it fails instead.
#'
#' @param dir         Directory to create the script in.
#' @param exit_status Exit status the script returns.
fake_aracne <- function(dir, exit_status = 0) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  network <- file.path(dir, "network.tsv")
  edges <- sim$planted_edges
  writeLines(c(paste("regulator.values", "target.values", "mi.values",
                     "scc.values", "count.values", "log.p.values", sep = "\t"),
               paste(edges$regulator, edges$target, 0.5, 0.5, 1, -10,
                     sep = "\t")),
             network)
  script <- file.path(dir, "fake_aracne3.sh")
  calls_log <- shQuote(file.path(dir, "calls.log"))
  writeLines(c(
    "#!/bin/bash",
    sprintf("printf '%%s\\n' \"$*\" >> %s", calls_log),
    sprintf("[ %d -eq 0 ] || exit %d", exit_status, exit_status),
    "while [ $# -gt 0 ]; do",
    "  if [ \"$1\" = \"-o\" ]; then out=\"$2\"; fi",
    "  shift",
    "done",
    sprintf("cp %s \"$out/consolidated-net_defaultid.tsv\"", shQuote(network))
  ), script)
  Sys.chmod(script, "755")
  script
}

#' Build a Runner on Fresh Output Directories
#'
#' @param bin  Path to the ARACNe3 executable to call.
#' @param root Directory under which all outputs are created.
#' @param ...  Passed to `Runner$new()` (e.g. `num_subnets`).
make_runner <- function(bin, root = tempfile("runner_"), mats = metacells$mats,
                        ...) {
  dirs <- file.path(root, c("out", "out/aracne_results", "out/viper_results"))
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  reg_file <- file.path(root, "tfs.txt")
  writeLines(sim$regulators, reg_file)
  Runner$new(mats, dirs[1], bin, dirs[2], list(tfs = reg_file),
             sim$seurat_obj, dirs[3], ...)
}

read_calls <- function(bin) {
  readLines(file.path(dirname(bin), "calls.log"))
}

test_that("ARACNe runs once per matrix and regulator set, on these matrices", {
  bin <- fake_aracne(tempfile("bin_"))
  runner <- make_runner(bin)
  runner$run_aracne()

  calls <- read_calls(bin)
  expect_length(calls, 3)
  for (mat_name in names(metacells$mats)) {
    input <- file.path(runner$aracne_output_path,
                       paste0(mat_name, "_for_aracne.tsv"))
    out_dir <- file.path(runner$aracne_output_path, paste0("tfs_", mat_name))
    expect_true(any(grepl(input, calls, fixed = TRUE)))
    expect_true(file.exists(file.path(out_dir,
                                      "consolidated-net_defaultid.tsv")))
  }
})

test_that("the ARACNe input file holds exactly the metacell matrix", {
  bin <- fake_aracne(tempfile("bin_"))
  runner <- make_runner(bin)
  runner$run_aracne()

  for (mat_name in names(metacells$mats)) {
    path <- file.path(runner$aracne_output_path,
                      paste0(mat_name, "_for_aracne.tsv"))
    header <- readLines(path, n = 1)
    mat <- metacells$mats[[mat_name]]
    # ARACNe3 counts samples from the tabs in the header line.
    expect_equal(lengths(regmatches(header, gregexpr("\t", header))),
                 ncol(mat))
    expect_true(startsWith(header, "\t"))
    expect_equal(read_metacell_tsv(path), mat, tolerance = 1e-12)
  }
})

test_that("threads, seed and num_subnets are passed to ARACNe3", {
  bin <- fake_aracne(tempfile("bin_"))
  make_runner(bin, threads = 2, seed = 7, num_subnets = 10)$run_aracne()
  calls <- read_calls(bin)
  expect_true(all(grepl("--threads 2 --seed 7 -x 10$", calls)))
})

test_that("num_subnets defaults to 1, ARACNe3's own default", {
  bin <- fake_aracne(tempfile("bin_"))
  runner <- make_runner(bin)
  runner$run_aracne()
  expect_equal(runner$num_subnets, 1)
  expect_true(all(grepl("-x 1$", read_calls(bin))))
})

test_that("files left in the output folder by older runs are not used", {
  bin <- fake_aracne(tempfile("bin_"))
  runner <- make_runner(bin)
  # What the previous version of the pipeline read and wrote.
  stale <- file.path(runner$base_output_path,
                     "metacell_clust-9-metaCells_all_all.txt.tsv")
  writeLines("stale", stale)
  stale_run <- file.path(runner$aracne_output_path, "tfs_metacell_clust-1")
  dir.create(file.path(stale_run, "subnets"), recursive = TRUE)
  writeLines("stale", file.path(stale_run, "subnets", "subnet9_defaultid.tsv"))

  runner$run_aracne()

  calls <- read_calls(bin)
  expect_length(calls, 3)
  expect_false(any(grepl("all_all", calls)))
  # The run directory is cleared before ARACNe writes into it.
  expect_false(file.exists(file.path(stale_run, "subnets",
                                     "subnet9_defaultid.tsv")))
})

test_that("paths containing spaces are passed safely", {
  root <- file.path(tempfile("runner_"), "dir with spaces")
  bin <- fake_aracne(file.path(root, "bin dir"))
  make_runner(bin, root = root)$run_aracne()
  expect_length(read_calls(bin), 3)
})

test_that("an ARACNe3 failure stops the run", {
  bin <- fake_aracne(tempfile("bin_"), exit_status = 3)
  expect_error(make_runner(bin)$run_aracne(), "ARACNe3 failed \\(exit status 3")
})

test_that("a missing ARACNe3 binary stops the run", {
  expect_error(make_runner(file.path(tempdir(), "no-such-aracne"))$run_aracne(),
               "ARACNe3 failed")
})

test_that("invalid inputs are rejected with clear errors", {
  bin <- fake_aracne(tempfile("bin_"))
  expect_error(make_runner(bin, num_subnets = 0), "num_subnets")
  expect_error(make_runner(bin, num_subnets = 2.5), "num_subnets")
  expect_error(make_runner(bin, threads = 0), "threads")
  expect_error(make_runner(bin, mats = list())$run_aracne(), "empty")
  expect_error(make_runner(bin, mats = unname(metacells$mats))$run_aracne(),
               "unique name")
  duplicated_genes <- metacells$mats[1]
  rownames(duplicated_genes[[1]])[2] <- rownames(duplicated_genes[[1]])[1]
  expect_error(make_runner(bin, mats = duplicated_genes)$run_aracne(),
               "duplicate gene names")
})

test_that("run_viper uses only this run's networks", {
  skip_if_not_installed("viper")
  bin <- fake_aracne(tempfile("bin_"))
  runner <- make_runner(bin)
  runner$run_aracne()
  # A leftover network from an older run, which would break VIPER if read.
  old_run <- file.path(runner$aracne_output_path,
                       "tfs_metacell_clust-1-metaCells")
  dir.create(old_run)
  writeLines("not a network", file.path(old_run,
                                        "consolidated-net_defaultid.tsv"))

  suppressMessages(suppressWarnings(capture.output(runner$run_viper())))

  regulon_files <- list.files(runner$base_output_path, "_pruned\\.rds$")
  expect_setequal(regulon_files,
                  paste0("tfs_", names(metacells$mats), "_pruned.rds"))
  viper_res <- readRDS(file.path(runner$viper_output_path,
                                 "viper_results.rds"))
  expect_setequal(rownames(viper_res), sim$regulators)
  expect_equal(colnames(viper_res), colnames(sim$seurat_obj))
})

test_that("run_viper explains when ARACNe has not been run", {
  skip_if_not_installed("viper")
  runner <- make_runner(fake_aracne(tempfile("bin_")))
  expect_error(runner$run_viper(), "run run_aracne\\(\\) first")
})
