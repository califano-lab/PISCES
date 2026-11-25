# Run the PISCES test suite from the repository root:
#
#   Rscript tests/run_tests.r
#
# The end-to-end ARACNe3 tests also run when ARACNE3_BIN points to an ARACNe3
# executable:
#
#   ARACNE3_BIN=/path/to/ARACNe3_app_release Rscript tests/run_tests.r
#
# Requires: testthat, R6, Seurat; viper for the VIPER tests (skipped without).

if (!file.exists(file.path("classes", "runner.r"))) {
  stop("Run this script from the PISCES repository root.")
}

testthat::test_dir(file.path("tests", "testthat"), reporter = "progress",
                   stop_on_failure = TRUE)
