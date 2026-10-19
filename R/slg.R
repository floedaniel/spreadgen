# Model SLG (paper: model C): simple logistic growth in every suitable cell, spatially implicit.

#' Logistic density in percent of K
#'
#' `N_t = 100 * N0 * exp(r t) / (100 + N0 * (exp(r t) - 1))`, the logistic used by all Robinet et al. (2012) models,
#' evaluated as `N0 / (N0 / 100 + (1 - N0 / 100) * exp(-r t))` (same value, no overflow for large `r t`).
#'
#' @param n0 Initial density, percent of K (0-100).
#' @param r Relative growth rate (yr-1).
#' @param t Time (years).
#' @return Density at `t`, percent of K.
#' @export
logistic_pct <- function(n0, r, t) n0 / (n0 / 100 + (1 - n0 / 100) * exp(-r * t))

#' Growth rate from observed invaded cells (model A / Rand-RR parameter)
#'
#' Inverts the logistic for the relative rate of increase of invaded cells:
#' `r = (1/t) ln( N_t / (N_max - N_t) * (N_max - N_0) / N_0 )`.
#'
#' @param n_t Invaded cells at time `t`.
#' @param n_max Number of suitable cells.
#' @param n_0 Invaded cells at t = 0.
#' @param t Years between the two counts.
#' @return r (yr-1). Example (D. virgifera): `r_from_invaded(358, 3104, 1, 18)` = 0.33.
#' @export
r_from_invaded <- function(n_t, n_max, n_0, t) {
  stopifnot(all(n_t > 0), all(n_t < n_max), all(n_0 > 0), all(n_0 < n_max), all(t > 0))
  log(n_t / (n_max - n_t) * (n_max - n_0) / n_0) / t
}

#' Simple logistic growth (SLG) model
#'
#' Every suitable cell starts at density `n0` (percent of K) and grows logistically with
#' `r_cell = GI * ln(lmax) / 100`; unsuitable cells stay 0. No dispersal.
#'
#' @param pra An `sg_pra` object ([prep_suitability()] or [prep_climex()]).
#' @param n0 Initial density, percent of K; a number or a SpatRaster on the PRA grid.
#'   For a founder population use [n0_from()] (e.g. 100 founders / K).
#' @param lmax Maximum yearly multiplication factor (where GI = 100); > 1.
#' @param t Years (vector) at which density is returned.
#' @param reach Density threshold (percent of K) for "reached" (default 1e-6).
#' @param abundant Density threshold (percent of K) for "abundant" (default 1; the paper's invaded threshold).
#' @return An object of class `sg_slg`: `density` (SpatRaster, one layer per `t`, percent of K), `summary`
#'   (data frame per `t`), `params`.
#' @export
spread_slg <- function(pra, n0, lmax, t = c(5, 10, 20, 50), reach = 1e-6, abundant = 1) {
  stopifnot(inherits(pra, "sg_pra"), lmax > 1, all(t >= 0), reach > 0, abundant >= reach)
  r <- pra$rast
  growth <- r[["gi"]] * log(lmax) / 100
  dens <- terra::rast(lapply(t, function(tt) logistic_pct(n0, growth, tt) * r[["ei"]]))
  names(dens) <- paste0("t", t)
  out <- list(density = dens, summary = density_summary(dens, r, t, reach, abundant),
              params = list(model = "SLG", n0 = if (is.numeric(n0)) n0 else "raster", lmax = lmax,
                            reach = reach, abundant = abundant))
  class(out) <- "sg_slg"
  out
}

#' Years until a cell reaches a density level under SLG
#'
#' Closed-form inverse of the logistic: `t = ln( N (100 - N0) / (N0 (100 - N)) ) / r`. Cells that never reach
#' `level` (unsuitable, GI = 0) are `Inf`; cells already at or above `level` are 0.
#'
#' @inheritParams spread_slg
#' @param level Density (percent of K, 0 < level < 100).
#' @return SpatRaster of years (NA where the PRA has no data).
#' @export
slg_time_to <- function(pra, n0, lmax, level = 1) {
  stopifnot(inherits(pra, "sg_pra"), lmax > 1, level > 0, level < 100)
  r <- pra$rast
  growth <- r[["gi"]] * log(lmax) / 100
  tt <- log(level * (100 - n0) / (n0 * (100 - level))) / growth
  tt <- terra::ifel(r[["ei"]] == 1 & growth > 0, tt, Inf)
  tt <- if (is.numeric(n0)) { if (n0 >= level) tt * 0 else tt } else terra::ifel(n0 >= level & r[["ei"]] == 1, 0, tt)
  names(tt) <- paste0("years_to_", level)
  terra::mask(tt, r[["ei"]])
}

#' Original SLG class counts (legacy, for regression against Prog_SpreadModule_V21-web.R)
#'
#' Reproduces `sum0, sum25, sum50, sum75, sum100` of the original `slg()`, including its boundary handling
#' (class 1: `0 < N <= 25`, class 2: `25 <= N < 50`, class 3: `50 <= N < 75`, class 4: `N >= 75`).
#' @param x An `sg_slg` object.
#' @return Data frame per `t`.
#' @export
slg_classes_v21 <- function(x) {
  v <- terra::values(x$density, mat = TRUE)
  v <- v[!is.na(v[, 1]), , drop = FALSE]
  data.frame(t = as.numeric(sub("^t", "", colnames(v))),
             sum0 = colSums(v == 0), sum25 = colSums(v > 0 & v <= 25), sum50 = colSums(v >= 25 & v < 50),
             sum75 = colSums(v >= 50 & v < 75), sum100 = colSums(v >= 75), row.names = NULL)
}

#' @export
print.sg_slg <- function(x, ...) {
  cat(sprintf("<sg_slg> n0 = %s %%, lmax = %g; thresholds reach %g %%, abundant %g %% of K\n",
              format(x$params$n0), x$params$lmax, x$params$reach, x$params$abundant))
  print(x$summary, row.names = FALSE, digits = 4)
  invisible(x)
}

# ---- internal ----

# Per-layer summary of a density stack (percent of K) over the PRA
density_summary <- function(dens, r, t, reach, abundant) {
  v <- terra::values(dens, mat = TRUE)
  ok <- !is.na(terra::values(r[["ei"]], mat = FALSE))
  suit <- ok & terra::values(r[["ei"]], mat = FALSE) == 1
  a <- terra::values(r[["area_km2"]], mat = FALSE)
  n_suit <- sum(suit)
  a_suit <- sum(a[suit])
  do.call(rbind, lapply(seq_along(t), function(i) {
    d <- v[, i]
    in_r <- ok & d > reach
    in_a <- ok & d > abundant
    data.frame(t = t[i], n_suitable = n_suit, area_suitable_km2 = a_suit,
               n_reach = sum(in_r), area_reach_km2 = sum(a[in_r]), pct_suitable_reach = 100 * sum(a[in_r]) / a_suit,
               n_abundant = sum(in_a), area_abundant_km2 = sum(a[in_a]),
               pct_suitable_abundant = 100 * sum(a[in_a]) / a_suit,
               mean_density_suitable = mean(d[suit]))
  }))
}
