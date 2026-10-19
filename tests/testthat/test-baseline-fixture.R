# The V21 baseline fixture (data-raw/run_v21_baseline.R) must reproduce the published D. virgifera results
# (Robinet et al. 2012, PLoS ONE 7: e43366, Fig. 6 areas; N_max = 3104 suitable cells).

b <- readRDS(test_path("fixtures", "v21_dvv_baseline.rds"))
key <- function(lat, lon) paste(round(lat, 3), round(lon, 3))
tab_area <- b$grid$area_km2[match(key(b$tab$lat, b$tab$lon), key(b$grid$lat, b$grid$lon))]

test_that("fixture has the expected structure", {
  expect_equal(nrow(b$tab), 6683)
  expect_equal(nrow(b$grid), 20995)
  expect_equal(dim(b$slg), c(6683, 18))
  expect_equal(dim(b$radial), c(6683, 18))
  expect_equal(dim(b$dispk), c(20995, 18))
  expect_false(anyNA(tab_area))
})

test_that("suitable cells match the paper (N_max = 3104)", {
  expect_equal(sum(b$tab$EI > 0), 3104)
  expect_equal(sum(b$grid$EI > 0), 3104)
})

test_that("published invaded areas are reproduced", {
  expect_equal(sum(tab_area[b$slg[, 8] > 1]), 59790, tolerance = 1e-5)         # model C, t = 8
  expect_equal(sum(tab_area[b$radial[, 18]]), 3578880, tolerance = 1e-6)       # model B, t = 18
  expect_equal(sum(b$grid$area_km2[b$dispk[, 18] > 1]), 2010382, tolerance = 1e-6)  # model D, t = 18
  expect_gt(mean(b$radial[b$tab$EI > 0, 18]), 0.5)                              # "more than half" invaded by 2010
})

test_that("baseline series are monotone in time", {
  expect_true(all(diff(colSums(b$radial)) >= 0))
  expect_true(all(apply(b$slg, 1, function(x) all(diff(x) >= 0))))
})
