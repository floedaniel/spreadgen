# Parallel runs use the INSTALLED package in PSOCK workers. Install first (from the folder above the package):
#   R CMD INSTALL spreadgen
# The tests skip when spreadgen is not installed or the installed version differs from the source.

installed_matches <- function() {
  lib <- suppressWarnings(tryCatch(utils::packageDescription("spreadgen", lib.loc = .libPaths()), error = function(e) NULL))
  is.list(lib) && identical(lib$Version, as.character(utils::packageVersion("spreadgen")))
}

# One unit of work: builds its own inputs (no terra objects cross processes), runs a small Monte Carlo
toy_species <- function(key, n_draws = 4) {
  r <- terra::rast(nrows = 15, ncols = 20, xmin = 4e6, xmax = 4.02e6, ymin = 3e6, ymax = 3.015e6, crs = "EPSG:3035")
  terra::values(r) <- rep(seq(0.2, 1, length.out = 20), times = 15)
  g <- spreadgen::as_grid(spreadgen::prep_suitability(r, 0.3))
  m <- spreadgen::mc_spread(g, "DK", list(RR = spreadgen::prior_pert(2, 3, 5, log = TRUE),
                                          p = spreadgen::prior_pert(3, 5, 9), lmax = spreadgen::prior_pert(3, 5, 8, log = TRUE),
                                          init = spreadgen::prior_pert(1, 10, 50, log = TRUE)),
                            n_draws, t = c(2, 4), entry = list(cells = 150L), seed = 1, key = key,
                            dk_args = list(cutoff = "none"))
  data.frame(key = key, n_ok = m$n_ok, median_cells = stats::median(m$metrics$n_cells[m$metrics$level == 1 & m$metrics$t == 4]),
             mean_prob = mean(m$prob[, 2], na.rm = TRUE), first_RR = m$draws$RR[1])
}

test_that("workers = 1 and workers = 2 give identical results (installed package, PSOCK)", {
  skip_on_cran()
  skip_if_not(installed_matches(), "spreadgen not installed, or installed version differs from the source")
  keys <- c("AAA", "BBB", "CCC")
  one <- run_batches(keys, toy_species, workers = 1, log_fun = function(...) NULL)
  two <- run_batches(keys, toy_species, workers = 2, log_fun = function(...) NULL)
  cols <- c("key", "n_ok", "median_cells", "mean_prob", "first_RR")
  expect_equal(one[, cols], two[, cols])
  expect_equal(length(unique(two$worker_pid)), 2)
  expect_false(Sys.getpid() %in% two$worker_pid)
})

test_that("parallel smoke test: 2 workers, 4 draws; failures become FAILED rows; on_batch sees every batch", {
  skip_on_cran()
  skip_if_not(installed_matches(), "spreadgen not installed, or installed version differs from the source")
  seen <- list()
  f <- function(key) if (key == "BAD") stop("boom") else toy_species(key)
  res <- run_batches(c("AAA", "BAD", "CCC"), f, workers = 2, on_batch = function(rows) seen[[length(seen) + 1]] <<- rows,
                     log_fun = function(...) NULL)
  expect_equal(nrow(res), 3)
  expect_equal(res$status[res$key == "BAD"], "FAILED: boom")
  expect_true(all(res$n_ok[res$key != "BAD"] == 4))
  expect_equal(length(seen), 2)                               # batches of 2 workers: (AAA, BAD), (CCC)
  expect_error(run_batches("AAA", toy_species, workers = 2, expected_version = "0.0.0", log_fun = function(...) NULL),
               "differs")
})
