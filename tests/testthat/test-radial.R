jaccard <- function(a, b) sum(a & b) / sum(a | b)

test_that("dist_v21: sphere of radius 20000/pi, rounding and clamping", {
  expect_equal(dist_v21(0, 0, 1, 0, round_digits = NULL), 20000 / 180, tolerance = 1e-9)
  expect_equal(dist_v21(0, 0, 1, 0), 111.109, tolerance = 1e-5)   # 7-digit cosine rounding: -1.8 km at 1 degree
  expect_equal(dist_v21(10, 50, 10, 50), 0)
  expect_equal(dist_v21(10, 50, 10, 50, round_digits = NULL), 0)
  expect_false(is.nan(dist_v21(10, 50, 10 + 1e-9, 50, round_digits = NULL)))
  expect_equal(dist_v21(0, 0, 180, 0), 20000)
})

test_that("RR reproduces the V21 baseline exactly (distance = 'v21', Jaccard = 1, t = 1..18)", {
  b <- v21()
  p <- pra_v21(b)
  x <- spread_radial(p, cbind(20.30, 44.82), RR = 80, t = 1:18, distance = "v21")
  cells <- terra::cellFromXY(p$rast, cbind(b$tab$lon, b$tab$lat))
  inv <- terra::values(x$invaded, mat = TRUE)[cells, ] == 1
  for (i in 1:18) expect_equal(jaccard(inv[, i], b$radial[, i]), 1, label = paste("Jaccard t =", i))
  expect_equal(x$summary$n_invaded, b$summary$radial$ninv)
  expect_equal(x$summary$pniche, b$summary$radial$pniche, tolerance = 1e-12)
  # published area at t = 18 (Fig. 6, model B), v21 cell areas
  expect_equal(x$summary$area_invaded_km2[18], 3578880, tolerance = 1e-4)
})

test_that("RR default (sphere, no rounding) and grid modes are close to V21", {
  b <- v21()
  p <- pra_v21(b)
  cells <- terra::cellFromXY(p$rast, cbind(b$tab$lon, b$tab$lat))
  for (m in c("sphere", "grid")) {
    x <- spread_radial(p, cbind(20.30, 44.82), RR = 80, t = c(5, 10, 18), distance = m)
    inv <- terra::values(x$invaded, mat = TRUE)[cells, ] == 1
    j <- sapply(1:3, function(i) jaccard(inv[, i], b$radial[, c(5, 10, 18)[i]]))
    expect_gt(min(j), if (m == "sphere") 0.999 else 0.97)
  }
})

test_that("RR: invaded set = suitable cells inside the disc; monotone; years_to_arrival consistent", {
  p <- prep_suitability(synth_suit(), 0.3)
  e <- cbind(4015500, 3010500)
  x <- spread_radial(p, e, RR = 2, t = 0:8, distance = "grid")
  v <- terra::values(x$invaded, mat = TRUE)
  ei <- terra::values(p$rast[["ei"]], mat = FALSE)
  xy <- terra::xyFromCell(p$rast, seq_len(terra::ncell(p$rast)))
  d <- sqrt((xy[, 1] - e[1])^2 + (xy[, 2] - e[2])^2) / 1000
  for (i in seq_along(0:8)) expect_equal(v[, i] == 1, ei == 1 & d <= 2 * (0:8)[i], label = paste("t =", (0:8)[i]))
  expect_true(all(apply(v, 1, function(z) all(diff(z) >= 0))))
  yr <- terra::values(x$years_to_arrival, mat = FALSE)
  expect_equal(yr[ei == 1], ceiling(d[ei == 1] / 2))
  expect_true(all(is.infinite(yr[ei == 0])))
  for (i in seq_along(0:8)) expect_equal(v[, i] == 1, yr <= (0:8)[i])
  expect_equal(sum(v[, 1]), 1)                         # t = 0: the entry cell only
})

test_that("RR: union of discs around entry cells only; unsuitable entry snaps", {
  p <- prep_suitability(synth_suit(), 0.3)
  one <- spread_radial(p, cbind(4015500, 3010500), RR = 1, t = 3, distance = "grid")
  two <- spread_radial(p, rbind(c(4015500, 3010500), c(4025500, 3005500)), RR = 1, t = 3, distance = "grid")
  oth <- spread_radial(p, cbind(4025500, 3005500), RR = 1, t = 3, distance = "grid")
  vv <- function(z) terra::values(z$invaded, mat = FALSE) == 1
  expect_equal(vv(two), vv(one) | vv(oth))
  # entry in an unsuitable cell (column 2) snaps to the first suitable column (10)
  s <- spread_radial(p, cbind(4001500, 3010500), RR = 1, t = 0, distance = "grid")
  expect_equal(terra::colFromCell(p$rast, s$entry), 10)
  expect_error(spread_radial(p, 1L, RR = 1))           # cell 1 is unsuitable
})

test_that("RR reach = 'land' goes around no-data barriers", {
  s <- synth_suit() * 0 + 1
  s[, 15] <- NA                                        # wall at column 15 ...
  s[20, 15] <- 1                                       # ... with a gap in the last row
  p <- prep_suitability(s, 0.5)
  e <- cbind(4010500, 3015500)                         # column 11, row 5
  disc <- spread_radial(p, e, RR = 1, t = 30, distance = "grid")
  land <- spread_radial(p, e, RR = 1, t = 30, reach = "land")
  target <- terra::cellFromRowCol(p$rast, 5, 19)
  d_disc <- terra::values(disc$distance_km, mat = FALSE)[target]
  d_land <- terra::values(land$distance_km, mat = FALSE)[target]
  expect_equal(d_disc, 8)
  expect_gt(d_land, 25)
  expect_lt(terra::values(disc$years_to_arrival, mat = FALSE)[target], terra::values(land$years_to_arrival, mat = FALSE)[target])
  expect_error(spread_radial(p, e, RR = 1, reach = "land", distance = "v21"))
})

test_that("Rand-RR: logistic count inside the radial range, nested in time, reproducible", {
  p <- prep_suitability(synth_suit(), 0.3)
  e <- cbind(4015500, 3010500)
  n_suit <- sum(terra::values(p$rast[["ei"]]) == 1, na.rm = TRUE)
  set.seed(3)
  x <- spread_radialrand(p, e, RR = 3, n0 = n0_from(1, n_suit), r = 0.6, t = c(2, 4, 6, 10), distance = "grid")
  rad <- spread_radial(p, e, RR = 3, t = c(2, 4, 6, 10), distance = "grid")
  expect_equal(x$summary$n_invaded,
               pmin(round(logistic_pct(n0_from(1, n_suit), 0.6, c(2, 4, 6, 10)) * n_suit / 100), rad$summary$n_invaded))
  v <- terra::values(x$invaded, mat = TRUE)
  vr <- terra::values(rad$invaded, mat = TRUE)
  expect_true(all(v[v == 1 & !is.na(v)] == vr[v == 1 & !is.na(v)]))      # inside the range
  expect_true(all(apply(v, 1, function(z) all(diff(z) >= 0)), na.rm = TRUE))
  set.seed(3)
  y <- spread_radialrand(p, e, RR = 3, n0 = n0_from(1, n_suit), r = 0.6, t = c(2, 4, 6, 10), distance = "grid")
  expect_equal(terra::values(y$invaded), terra::values(x$invaded))
})
