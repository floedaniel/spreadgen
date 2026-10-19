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

#' SLG model core on a plain-data grid
#'
#' @param g An `sg_grid` object ([as_grid()]).
#' @param n0 Initial density (percent of K): a number or a vector over all cells.
#' @param lmax Maximum yearly multiplication factor (> 1).
#' @param t Years (vector).
#' @return Matrix cells x `t` of density (percent of K; 0 in unsuitable cells, NA where the grid has no data).
#' @export
core_slg <- function(g, n0, lmax, t) {
  stopifnot(inherits(g, "sg_grid"), lmax > 1, all(t >= 0))
  growth <- g$gi * log(lmax) / 100
  m <- vapply(t, function(tt) logistic_pct(n0, growth, tt) * g$ei, numeric(length(g$ei)))
  matrix(m, ncol = length(t))
}

#' Years until each cell reaches a density level under SLG (core)
#'
#' Closed-form inverse of the logistic: `t = ln( N (100 - N0) / (N0 (100 - N)) ) / r`; `Inf` where the cell never
#' reaches `level` (unsuitable, GI = 0), 0 where it starts at or above it.
#' @inheritParams core_slg
#' @param level Density (percent of K, 0 < level < 100).
#' @return Numeric vector over all cells (NA where the grid has no data).
#' @export
core_slg_years <- function(g, n0, lmax, level) {
  stopifnot(lmax > 1, level > 0, level < 100)
  growth <- g$gi * log(lmax) / 100
  n0 <- rep_len(n0, length(g$ei))
  y <- log(level * (100 - n0) / (n0 * (100 - level))) / growth
  y[!(g$ei == 1 & growth > 0)] <- Inf
  y[g$ei == 1 & n0 >= level] <- 0
  y[is.na(g$ei)] <- NA
  y
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
#' @param levels Density levels (percent of K, strict `>`) for the summary; default `c(0.1, 1, 10)`, headline 1.
#' @return An object of class `sg_slg`: `density` (SpatRaster, one layer per `t`, percent of K), `summary`
#'   (data frame per `t` x level), `params`.
#' @export
spread_slg <- function(pra, n0, lmax, t = c(5, 10, 20, 50), levels = c(0.1, 1, 10)) {
  stopifnot(inherits(pra, "sg_pra"), all(levels > 0), all(levels < 100))
  g <- as_grid(pra)
  n0v <- if (inherits(n0, "SpatRaster")) terra::values(n0, mat = FALSE) else n0
  dens <- core_slg(g, n0v, lmax, t)
  out <- list(density = grid_to_rast(g, dens, paste0("t", t)), summary = level_summary(dens, g, t, levels),
              params = list(model = "SLG", n0 = if (is.numeric(n0)) n0 else "raster", lmax = lmax, levels = levels))
  class(out) <- "sg_slg"
  out
}

#' Years until a cell reaches a density level under SLG
#'
#' SpatRaster wrapper of [core_slg_years()].
#' @inheritParams spread_slg
#' @param level Density (percent of K, 0 < level < 100).
#' @return SpatRaster of years (`Inf` = never; NA where the PRA has no data).
#' @export
slg_time_to <- function(pra, n0, lmax, level = 1) {
  stopifnot(inherits(pra, "sg_pra"))
  g <- as_grid(pra)
  n0v <- if (inherits(n0, "SpatRaster")) terra::values(n0, mat = FALSE) else n0
  grid_to_rast(g, core_slg_years(g, n0v, lmax, level), paste0("years_to_", level))
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
  cat(sprintf("<sg_slg> n0 = %s %%, lmax = %g; levels %s %% of K\n", format(x$params$n0), x$params$lmax,
              paste(x$params$levels, collapse = " / ")))
  print(x$summary, row.names = FALSE, digits = 4)
  invisible(x)
}
