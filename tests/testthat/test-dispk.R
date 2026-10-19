# ---- kernel ----

test_that("2Dt kernel integrates to 1 over the plane; tail formula matches", {
  for (p in c(1.5, 2, 5, 50)) {
    f <- function(r) 2 * pi * r * kernel_2dt(r, u = 10, p = p)
    expect_equal(integrate(f, 0, Inf, rel.tol = 1e-10)$value, 1, tolerance = 1e-6, label = paste("p =", p))
    expect_equal(integrate(f, 30, Inf, rel.tol = 1e-10)$value, kernel_tail_2dt(30, 10, p), tolerance = 1e-6)
  }
})

test_that("2Dt kernel: p large -> normal with SD u; p -> 1+ -> Cauchy-like tail; no overflow", {
  u <- 7
  r <- c(0, 3, 7, 14, 21)
  expect_equal(kernel_2dt(r, u, 1e5), exp(-r^2 / (2 * u^2)) / (2 * pi * u^2), tolerance = 1e-3)
  expect_true(is.finite(kernel_2dt(10, 5, 1000)))      # gamma() would overflow here
  expect_false(is.finite(gamma(1001 / 2)))
  # tail exponent -> r^-(p+1): local log-slope at large r
  slope <- function(p) diff(log(kernel_2dt(c(1e4, 2e4), 1, p))) / log(2)
  expect_equal(slope(1.01), -2.01, tolerance = 1e-3)
  expect_gt(kernel_tail_2dt(100, 1, 1.01), 0.9)        # almost all mass beyond 100 u
  expect_lt(kernel_tail_2dt(100, 1, 50), 1e-50)
  expect_error(kernel_2dt(1, 1, 1))
})

test_that("sphere correction is 1 at r = 0 and grows with distance", {
  expect_equal(kernel_2dt(0, 50, 5, sphere = TRUE), kernel_2dt(0, 50, 5))
  expect_gt(kernel_2dt(3000, 50, 5, sphere = TRUE) / kernel_2dt(3000, 50, 5), 1.03)
})

# ---- reference vs V21 ----

dvv_dispk_ref <- function(t) {
  b <- v21()
  p <- pra_v21(b)
  pres <- b$info$params$presence
  init <- init_density(p, as.matrix(pres[, 1:2]), pres[, 3])
  x <- spread_dispk_reference(p, init, lmax = 40, u = 80, p = 5, t = t, distance = "v21")
  key <- function(lat, lon) paste(round(lat, 3), round(lon, 3))
  xy <- terra::xyFromCell(p$rast, seq_len(terra::ncell(p$rast)))
  idx <- match(key(b$grid$lat, b$grid$lon), key(xy[, 2], xy[, 1]))
  d <- terra::values(x$density, mat = TRUE)[idx, , drop = FALSE]
  d[is.na(d)] <- 0                                          # V21 grid cells without CLIMEX data (sea) are 0
  list(x = x, d = d, ref = b$dispk[, t, drop = FALSE])
}

test_that("DK reference reproduces V21 densities (t = 1..2, lon/lat, v21 areas)", {
  z <- dvv_dispk_ref(1:2)
  expect_equal(unname(z$d), unname(z$ref), tolerance = 1e-10)
})

test_that("DK reference reproduces V21 for t = 1..18 (slow; SPREADGEN_SLOW_TESTS=true)", {
  skip_if_not(identical(Sys.getenv("SPREADGEN_SLOW_TESTS"), "true"), "slow test")
  z <- dvv_dispk_ref(1:18)
  expect_equal(unname(z$d), unname(z$ref), tolerance = 1e-10)
  s1 <- z$x$summary[z$x$summary$level == 1 & z$x$summary$t == 18, ]
  expect_equal(s1$n_cells, 976)
  expect_equal(s1$area_km2, 2010382, tolerance = 1e-6)
})

# ---- FFT vs reference (same projected grid) ----

synth_pra <- function(nr = 20, nc = 30, holes = TRUE) {
  s <- synth_suit(nr, nc)
  s <- s * 0.6 + 0.4
  if (holes) {
    s[3:6, 20:23] <- 0.05                                   # unsuitable block
    s[15, 5:10] <- NA                                       # no data
  }
  prep_suitability(s, 0.1)
}

test_that("FFT equals the reference implementation to 1e-8 (no cutoff)", {
  p <- synth_pra()
  init <- init_density(p, cbind(4005500, 3010500), 0.01)
  ref <- spread_dispk_reference(p, init, lmax = 6, u = 3, p = 4, t = 1:6)
  fft <- spread_dispk(p, init, lmax = 6, u = 3, p = 4, t = 1:6, cutoff = "none")
  expect_equal(terra::values(fft$density), terra::values(ref$density), tolerance = 1e-8)
  expect_equal(terra::values(fft$years_to), terra::values(ref$years_to))
  expect_equal(fft$mass$lost_unsuitable, ref$mass$lost_unsuitable, tolerance = 1e-8)
  expect_equal(fft$mass$lost_domain + fft$mass$lost_kernel, ref$mass$lost_domain, tolerance = 1e-8)
  # P < 1 (stayers) agrees as well
  ref2 <- spread_dispk_reference(p, init, lmax = 6, u = 3, p = 4, t = 1:3, P = 0.7)
  fft2 <- spread_dispk(p, init, lmax = 6, u = 3, p = 4, t = 1:3, P = 0.7, cutoff = "none")
  expect_equal(terra::values(fft2$density), terra::values(ref2$density), tolerance = 1e-8)
})

test_that("mass balance: no growth, all suitable, no cutoff", {
  p <- prep_suitability(synth_suit(40, 40) * 0 + 1, 0.5)
  init <- init_density(p, cbind(4020500, 3020500), 50)
  x <- spread_dispk(p, init, lmax = 1, u = 2, p = 5, t = 1:4, cutoff = "none")
  m <- x$mass
  # identity: dispersing x kernel weight sum = kept + left the grid (+ unsuitable = 0)
  expect_equal(m$dispersing * x$params$kernel_weight_sum, m$mass_after_dispersal + m$lost_domain + m$lost_unsuitable,
               tolerance = 1e-10)
  expect_equal(m$mass_after_growth, m$mass_after_dispersal)
  expect_equal(m$lost_unsuitable, rep(0, 4))
  # closed domain in practice: kernel well inside the grid -> mass conserved before truncation
  expect_lt(m$lost_domain[1] / m$dispersing[1], kernel_tail_2dt(19, 2, 5))   # source 19.5 km from the edge
  expect_equal(x$params$kernel_weight_sum, 1, tolerance = 5e-3)
  expect_equal(m$mass_after_dispersal[1], 50, tolerance = 5e-3)
})

test_that("auto cutoff meets its criterion; time-to-level insensitive to the cutoff", {
  p <- synth_pra(30, 40, holes = FALSE)
  init <- init_density(p, cbind(4002500, 3015500), 0.5)
  a <- spread_dispk(p, init, lmax = 5, u = 2, p = 3, t = 1:8, cutoff = "auto")
  n <- spread_dispk(p, init, lmax = 5, u = 2, p = 3, t = 1:8, cutoff = "none")
  R_c <- a$params$R_c_requested_km
  n_suit <- sum(terra::values(p$rast[["ei"]]) == 1, na.rm = TRUE)
  expect_equal(n_suit * 100 * kernel_2dt(R_c, 2, 3) * 1, 0.1 * min(a$params$levels), tolerance = 1e-8)
  expect_equal(names(a$years_to), c("years_to_0.1", "years_to_1", "years_to_10"))
  expect_equal(terra::values(a$years_to), terra::values(n$years_to))
  k3 <- spread_dispk(p, init, lmax = 5, u = 2, p = 3, t = 1:8, cutoff = "k", k = 3)
  expect_equal(k3$params$R_c_km, 6)
  expect_equal(k3$params$tail_mass_beyond_cutoff, kernel_tail_2dt(6, 2, 3))
  expect_gt(k3$mass$lost_kernel[1], 0)
})

test_that("auto cutoff with a thin-tailed kernel truncates, and the error stays below 0.1 x lowest level", {
  p <- synth_pra(40, 60, holes = FALSE)
  init <- init_density(p, cbind(4002500, 3020500), 50)
  a <- spread_dispk(p, init, lmax = 3, u = 2, p = 60, t = 1:6, cutoff = "auto")
  n <- spread_dispk(p, init, lmax = 3, u = 2, p = 60, t = 1:6, cutoff = "none")
  expect_false(a$params$cutoff_capped)
  expect_lt(a$params$R_c_km, 40)
  expect_lt(max(abs(terra::values(a$density) - terra::values(n$density)), na.rm = TRUE), 0.1 * min(a$params$levels))
  expect_equal(terra::values(a$years_to), terra::values(n$years_to))
})

test_that("skip below grid, weight-sum warning, levels checked, lon/lat refused", {
  p <- synth_pra()
  init <- init_density(p, cbind(4005500, 3010500), 1)
  expect_warning(s <- spread_dispk(p, init, lmax = 5, u = 1.5, p = 5, t = 1:2), "spread below grid")
  expect_equal(s$status, "spread below grid")
  expect_warning(spread_dispk(p, init, lmax = 5, u = 0.3, p = 5, t = 1, min_u_cells = 0, cutoff = "k"), "sum to")
  expect_error(spread_dispk(p, init, lmax = 5, u = 3, p = 5, t = 1, levels = c(0, 1)))
  expect_false("reach" %in% names(formals(spread_dispk)))
  expect_error(spread_dispk(pra_v21(), init, lmax = 5, u = 80, p = 5), "projected")
})

test_that("cell-integrated kernel has a better mass balance than point values at small u", {
  p <- synth_pra(40, 40, holes = FALSE)
  init <- init_density(p, cbind(4020500, 3020500), 1)
  pt <- suppressWarnings(spread_dispk(p, init, lmax = 1, u = 0.6, p = 5, t = 1, min_u_cells = 0, cutoff = "k", k = 20))
  ce <- spread_dispk(p, init, lmax = 1, u = 0.6, p = 5, t = 1, min_u_cells = 0, cutoff = "k", k = 20, kernel = "cell")
  expect_gt(abs(pt$params$kernel_weight_sum - 1), abs(ce$params$kernel_weight_sum - 1))
})

test_that("time_to_event: first cell and area shares, censored at the horizon", {
  p <- prep_suitability(synth_suit(2, 5) * 0 + 1, 0.5)          # 10 cells, 1 km2 each
  y <- terra::setValues(p$rast[["ei"]], c(0, 1, 2, 3, 4, 5, 6, Inf, Inf, Inf))
  te <- time_to_event(y, p, pct = c(10, 50, 80), horizon = 5)
  expect_equal(te$first_year, 0)
  expect_equal(te$year_pct_10, 0)
  expect_equal(te$year_pct_50, 4)
  expect_true(is.infinite(te$year_pct_80) && te$censored_year_pct_80)
  reg <- terra::setValues(p$rast[["ei"]], c(rep(NA, 5), rep(1, 5)))
  te2 <- time_to_event(y, p, target = reg, pct = 10, horizon = 5)
  expect_equal(te2$first_year, 5)
  expect_equal(te2$n_suitable_target, 5)
})
