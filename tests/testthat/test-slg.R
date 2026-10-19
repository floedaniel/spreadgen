test_that("logistic_pct is the closed-form logistic and saturates at 100", {
  expect_equal(logistic_pct(5, 0.5, 0), 5)
  n <- 1e-3; r <- 0.7; t <- 3
  expect_equal(logistic_pct(n, r, t), 100 / (1 + (100 / n - 1) * exp(-r * t)))
  expect_equal(logistic_pct(1e-7, log(40), 100), 100)
  expect_true(all(diff(logistic_pct(1e-5, 1, 0:30)) >= 0))
})

test_that("SLG on a single cell matches the closed form", {
  r1 <- terra::rast(nrows = 1, ncols = 1, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000, crs = "EPSG:3035", vals = 0.8)
  p <- prep_suitability(r1, 0.1)
  x <- spread_slg(p, n0 = 1e-4, lmax = 10, t = c(1, 5, 10))
  expect_equal(unname(terra::values(x$density)[1, ]), logistic_pct(1e-4, 80 * log(10) / 100, c(1, 5, 10)))
})

test_that("SLG is monotone in t, lmax and GI; unsuitable cells stay 0", {
  p <- prep_suitability(synth_suit(), 0.2)
  v <- function(...) terra::values(spread_slg(p, ...)$density, mat = TRUE)
  d <- v(n0 = 1e-3, lmax = 5, t = 0:20)
  ok <- !is.na(d[, 1])
  expect_true(all(apply(d[ok, ], 1, function(z) all(diff(z) >= 0))))
  expect_true(all(v(n0 = 1e-3, lmax = 10, t = 5) >= v(n0 = 1e-3, lmax = 5, t = 5)))
  row1 <- d[1:30, 6]                                       # one row, GI increasing west -> east
  expect_true(all(diff(row1[row1 > 0]) > 0))
  ei <- terra::values(p$rast[["ei"]], mat = FALSE)
  expect_true(all(d[ei == 0, ] == 0))
  expect_equal(max(v(n0 = 1e-3, lmax = 5, t = 500)), 100)
})

test_that("SLG reproduces the V21 baseline exactly (D. virgifera, t = 1..18)", {
  b <- v21()
  p <- pra_v21(b)
  x <- spread_slg(p, n0 = b$info$params$N0_slg, lmax = b$info$params$lmax, t = 1:18)
  cells <- terra::cellFromXY(p$rast, cbind(b$tab$lon, b$tab$lat))
  d <- terra::values(x$density, mat = TRUE)[cells, ]
  expect_equal(unname(d), unname(b$slg), tolerance = 1e-12)
  cl <- slg_classes_v21(x)
  expect_equal(cl[, -1], b$summary$slg[, c("sum0", "sum25", "sum50", "sum75", "sum100")], ignore_attr = TRUE)
  # published area at t = 8 (Fig. 6, model C): 59,790 km2 above 1 % of K
  expect_equal(x$summary$area_abundant_km2[8], 59790, tolerance = 1e-4)
  expect_equal(x$summary$n_suitable[1], 3104)
})

test_that("slg_time_to inverts spread_slg", {
  p <- prep_suitability(synth_suit(), 0.2)
  tt <- slg_time_to(p, n0 = 1e-4, lmax = 20, level = 1)
  v <- terra::values(tt, mat = FALSE)
  ei <- terra::values(p$rast[["ei"]], mat = FALSE)
  expect_true(all(is.infinite(v[ei == 0])))
  i <- which(ei == 1)[5]
  gi <- terra::values(p$rast[["gi"]], mat = FALSE)[i]
  expect_equal(logistic_pct(1e-4, gi * log(20) / 100, v[i]), 1, tolerance = 1e-10)
  expect_true(all(terra::values(slg_time_to(p, n0 = 5, lmax = 20, level = 1))[ei == 1] == 0))
})

test_that("summary thresholds: reach <= abundant counts", {
  p <- prep_suitability(synth_suit(), 0.2)
  s <- spread_slg(p, n0 = 1e-5, lmax = 8, t = c(2, 6, 12))$summary
  expect_true(all(s$n_reach >= s$n_abundant))
  expect_true(all(diff(s$pct_suitable_abundant) >= 0))
  expect_true(all(s$pct_suitable_reach <= 100))
  expect_error(spread_slg(p, n0 = 1e-5, lmax = 8, reach = 2, abundant = 1))
})
