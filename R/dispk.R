# Model DK (paper: model D): 2Dt dispersal kernel + logistic growth, deterministic integrodifference model.
# Per year: dispersal (kernel value at the centre distance x landing-cell area, no renormalisation), stayers,
# death in unsuitable cells, logistic growth with r = GI ln(lmax) / 100. Density in percent of K.

#' 2Dt dispersal kernel (Clark et al. 1999)
#'
#' `f(r) = Gamma((p+1)/2) / (Gamma((p-1)/2) u^2 p pi) * (1 + r^2 / (p u^2))^(-(p+1)/2)` (km-2), evaluated with
#' `lgamma` (no overflow for large `p`). p -> 1: Cauchy-like fat tail; p -> Inf: bivariate normal with SD `u`.
#'
#' @param r Distance (km).
#' @param u Scale (km), > 0.
#' @param p Shape (degrees of freedom; the paper's nu), > 1.
#' @param sphere Apply the original sphere correction `1/s`, `s = sin(r/rho) / (r/rho)`, rho = 20000/pi (lon/lat
#'   grids).
#' @return Probability density per km2.
#' @export
kernel_2dt <- function(r, u, p, sphere = FALSE) {
  stopifnot(length(u) == 1, u > 0, length(p) == 1, p > 1)
  f <- exp(lgamma((p + 1) / 2) - lgamma((p - 1) / 2)) / (u^2 * p * pi) * (1 + r^2 / (p * u^2))^(-(p + 1) / 2)
  if (sphere) {
    rho <- 20000 / pi
    s <- sin(r / rho) / (r / rho)
    s[r == 0] <- 1
    f <- f / s
  }
  f
}

#' Probability mass of the 2Dt kernel beyond a distance
#'
#' `P(R > R_c) = (1 + R_c^2 / (p u^2))^(-(p-1)/2)` (planar kernel).
#' @inheritParams kernel_2dt
#' @param R_c Cutoff distance (km).
#' @return Tail mass (0-1).
#' @export
kernel_tail_2dt <- function(R_c, u, p) (1 + R_c^2 / (p * u^2))^(-(p - 1) / 2)

#' Initial density raster from entry points
#'
#' @param pra An `sg_pra` object.
#' @param entry Coordinates (two columns) or PRA cell numbers.
#' @param density Density at each entry (percent of K); recycled. Several points in one cell: mean (as the
#'   original's 3-column presence file).
#' @param crs CRS of coordinates when it differs from the PRA grid.
#' @param to Snap to `"any"` cell with data (default, as the original presence file) or `"suitable"`.
#' @return SpatRaster (percent of K, 0 elsewhere, NA where the PRA has no data).
#' @export
init_density <- function(pra, entry, density, crs = NULL, to = c("any", "suitable")) {
  to <- match.arg(to)
  cells <- if (is.null(dim(entry))) as.integer(entry) else snap_entry(pra, entry, crs = crs, to = to)$cell
  density <- rep_len(density, length(cells))
  v <- ifelse(is.na(terra::values(pra$rast[["ei"]], mat = FALSE)), NA, 0)
  m <- tapply(density, cells, mean)
  v[as.integer(names(m))] <- m
  out <- terra::setValues(pra$rast[["ei"]], v)
  names(out) <- "density"
  out
}

#' Dispersal kernel model, literal reference implementation
#'
#' Port of `dispk()` in `Prog_SpreadModule_V21-web.R`: loops over source cells (suitable or initially present)
#' and adds `f(d) * area_j * P * pop_i` to every cell `j` with data. Slow (O(sources x cells) per year); used for
#' regression and to check [spread_dispk()]. Works on lon/lat and projected grids.
#'
#' @param pra An `sg_pra` object.
#' @param init Initial density, SpatRaster from [init_density()] (percent of K).
#' @param lmax Maximum yearly multiplication factor (> 1; 1 = no growth).
#' @param u,p Kernel scale (km) and shape, see [kernel_2dt()].
#' @param t Years at which densities are returned (simulation runs yearly to `max(t)`).
#' @param P Proportion dispersing (default 1, as the tutorial advises).
#' @param distance `"auto"` (lon/lat: `"sphere"`, projected: `"euclid"`), `"v21"` (great circle with the
#'   original 7-digit cosine rounding), `"sphere"` (great circle, no rounding) or `"euclid"` (planar).
#' @param sphere_correction Apply the kernel's `1/s` sphere correction; default TRUE on lon/lat grids.
#' @param levels Density levels (percent of K, strict `>`) for `years_to` and the summary; default
#'   `c(0.1, 1, 10)`, headline 1 (the paper's invaded threshold).
#' @return An object of class `sg_dispk` (see [spread_dispk()]).
#' @export
spread_dispk_reference <- function(pra, init, lmax, u, p, t = 1:10, P = 1, distance = c("auto", "v21", "sphere", "euclid"),
                                   sphere_correction = NULL, levels = c(0.1, 1, 10)) {
  distance <- match.arg(distance)
  g <- as_grid(pra)
  lonlat <- g$lonlat
  if (distance == "auto") distance <- if (lonlat) "sphere" else "euclid"
  if (distance != "euclid" && !lonlat) stop("distance '", distance, "' needs a lon/lat grid")
  if (is.null(sphere_correction)) sphere_correction <- lonlat
  check_kernel_params(lmax, p, P)
  check_levels(levels)
  ok <- which(!is.na(g$ei))
  xy <- grid_xy(g, ok)
  if (!lonlat) xy <- xy / 1000
  area <- g$area_km2[ok]
  suit <- g$ei[ok] == 1
  growth <- g$gi[ok] * log(lmax) / 100
  pop <- terra::values(init, mat = FALSE)[ok]
  pop[is.na(pop)] <- 0
  sources <- which(suit | pop > 0)
  digits <- if (distance == "v21") 7 else NULL
  t_max <- max(t)
  keep <- list()
  yrs <- matrix(Inf, length(ok), length(levels))
  for (j in seq_along(levels)) yrs[pop > levels[j], j] <- 0
  mass <- vector("list", t_max)
  for (k in seq_len(t_max)) {
    popdisp <- numeric(length(ok))
    for (i in sources[pop[sources] > 0]) {
      R <- if (distance == "euclid") sqrt((xy[, 1] - xy[i, 1])^2 + (xy[, 2] - xy[i, 2])^2)
           else dist_v21(xy[i, 1], xy[i, 2], xy[, 1], xy[, 2], digits)
      popdisp <- popdisp + kernel_2dt(R, u, p, sphere = sphere_correction) * area * P * pop[i]
    }
    m_start <- sum(pop)
    landed <- sum(popdisp)
    popdisp <- popdisp + (1 - P) * pop
    lost_unsuit <- sum(popdisp[!suit])
    popdisp[!suit] <- 0
    pop <- 100 * popdisp * exp(growth) / (100 + popdisp * (exp(growth) - 1))
    mass[[k]] <- data.frame(year = k, mass_start = m_start, dispersing = P * m_start, lost_kernel = NA_real_,
                            lost_domain = P * m_start - landed, lost_unsuitable = lost_unsuit,
                            mass_after_dispersal = sum(popdisp), mass_after_growth = sum(pop))
    for (j in seq_along(levels)) yrs[is.infinite(yrs[, j]) & pop > levels[j], j] <- k
    if (k %in% t) keep[[as.character(k)]] <- pop
  }
  full <- function(x) { v <- rep(NA_real_, length(g$ei)); v[ok] <- x; v }
  res <- list(density = sapply(as.character(t), function(kk) full(keep[[kk]])),
              years = apply(yrs, 2, full), mass = do.call(rbind, mass), status = "ok",
              params = list(model = "DK", impl = "reference", lmax = lmax, u = u, p = p, P = P, distance = distance,
                            sphere_correction = sphere_correction, horizon = t_max, levels = levels))
  dispk_result(res, g, t)
}

#' Dispersal kernel (DK) model core: FFT convolution on a plain-data equal-area grid
#'
#' One dispersal step is a 2-D convolution of the density with the kernel value at the centre distance,
#' `f(r)`, multiplied by the landing-cell area (as the original; no renormalisation; on an equal-area grid this is
#' `f(r) h^2`), computed by zero-padded FFT (no wrap-around). Mass dispersing into unsuitable / no-data cells,
#' beyond the grid, or beyond the kernel cutoff is lost and reported per year.
#'
#' The kernel cutoff `R_c`: `"auto"` (default) is the distance beyond which all suitable cells together, each fully
#' occupied (100 % of K), add less than `tail_safety * min(levels)` to a landing cell,
#' i.e. `n_suitable * 100 f(R_c) h^2 <= tail_safety * min(levels)` (conservative); `"k"` uses `R_c = k * u`;
#' `"none"` uses the whole grid (no truncation inside it).
#'
#' @param g An `sg_grid` object on a projected, equal-area grid with square cells (e.g. EPSG:3035).
#' @param init Initial density (percent of K), vector over all cells.
#' @inheritParams spread_dispk_reference
#' @param cutoff `"auto"`, `"k"` or `"none"`.
#' @param k Cutoff factor for `cutoff = "k"`.
#' @param tail_safety Factor for `cutoff = "auto"` (default 0.1).
#' @param kernel `"point"` (default, original discretisation) or `"cell"` (kernel averaged over `cell_sub` x
#'   `cell_sub` points per landing cell; deviation, better mass balance when `u` is near the cell size).
#' @param cell_sub Sub-points per side for `kernel = "cell"`.
#' @param min_u_cells Minimum `u` in cells; below it the model is skipped (`status = "spread below grid"`).
#' @param keep_density Return the density matrix (cells x `t`); FALSE saves memory in Monte Carlo runs.
#' @return List: `density` (matrix cells x `t`, percent of K, or NULL), `years` (matrix cells x levels: first year
#'   above each level, `Inf` = not within the horizon), `mass` (per year, percent-of-K x cells: `mass_start`,
#'   `dispersing`, `lost_kernel` (1 - kernel weight sum: truncation and discretisation), `lost_domain` (landed outside
#'   the grid or in no-data cells), `lost_unsuitable` (landed in unsuitable cells with data),
#'   `mass_after_dispersal`, `mass_after_growth`), `params`, `status` (`"ok"` or the skip reason).
#' @export
core_dispk <- function(g, init, lmax, u, p, t = c(5, 10, 20, 50), P = 1, levels = c(0.1, 1, 10),
                       cutoff = c("auto", "k", "none"), k = 10, tail_safety = 0.1, kernel = c("point", "cell"),
                       cell_sub = 5, min_u_cells = 2, keep_density = TRUE) {
  cutoff <- match.arg(cutoff)
  kernel <- match.arg(kernel)
  stopifnot(inherits(g, "sg_grid"))
  if (g$lonlat) stop("the FFT kernel model needs a projected equal-area grid; use spread_dispk_reference() on lon/lat")
  if (abs(g$res[1] - g$res[2]) > 1e-6 * g$res[1]) stop("the FFT kernel model needs square cells")
  h <- g$cell_km
  check_kernel_params(lmax, p, P)
  check_levels(levels)
  t_max <- max(t)
  base_params <- list(model = "DK", impl = "fft", lmax = lmax, u = u, p = p, P = P, cell_km = h, horizon = t_max,
                      levels = levels)
  if (u / h < min_u_cells) {
    return(list(density = NULL, years = NULL, mass = NULL, status = "spread below grid", params = base_params))
  }
  nr <- g$nrow
  nc <- g$ncol
  tomat <- function(v) matrix(v, nr, nc, byrow = TRUE)
  ei <- tomat(g$ei)
  has_data <- !is.na(ei)
  suit <- has_data & ei == 1
  # FFT round-off is ~1e-12 of the maximum density: the cutoff criterion never goes below 1e-8 % of K
  target <- tail_safety * max(min(levels), 1e-8) / max(1, sum(suit))
  R_c <- switch(cutoff, none = Inf, k = k * u, auto = auto_cutoff(u, p, h, target))
  # kernel half-sizes per axis: offsets beyond the grid cannot land inside it
  n_c <- if (is.finite(R_c)) ceiling(R_c / h) else Inf
  n_r <- min(n_c, nr - 1)
  n_x <- min(n_c, nc - 1)
  # kernel density values f (km-2); landing mass = (f * pop) convolved, times the landing-cell area (as the original)
  w <- kernel_matrix(u, p, h, n_r, n_x, if (cutoff == "none") Inf else R_c, kernel, cell_sub)
  w_sum <- sum(w) * h^2
  if (w_sum > 1.05) warning(sprintf("kernel weights sum to %.3f (> 1.05): u is small relative to the cell size; consider kernel = 'cell'", w_sum))
  # circular convolution is exact inside the grid when the padded size is >= grid + half-kernel
  Nr <- stats::nextn(nr + n_r)
  Nc <- stats::nextn(nc + n_x)
  Kf <- stats::fft(embed_kernel(w, n_r, n_x, Nr, Nc))
  R_used <- if (is.finite(R_c)) min(R_c, h * sqrt(n_r^2 + n_x^2)) else h * sqrt(n_r^2 + n_x^2)

  area <- tomat(g$area_km2)
  area[is.na(area)] <- 0
  growth <- tomat(g$gi) * log(lmax) / 100
  growth[is.na(growth)] <- 0
  pop <- tomat(init)
  pop[is.na(pop)] <- 0
  yrs <- lapply(levels, function(lv) ifelse(pop > lv, 0, Inf))
  keep <- list()
  mass <- vector("list", t_max)
  A <- matrix(0, Nr, Nc)
  t0 <- proc.time()[["elapsed"]]
  for (yr in seq_len(t_max)) {
    m_start <- sum(pop)
    A[] <- 0
    A[seq_len(nr), seq_len(nc)] <- P * pop
    full <- Re(stats::fft(stats::fft(A) * Kf, inverse = TRUE)) / (Nr * Nc)
    inside <- pmax(full[seq_len(nr), seq_len(nc)], 0) * area
    popdisp <- inside + (1 - P) * pop
    lost_unsuit <- sum(popdisp[has_data & !suit])
    popdisp[!suit] <- 0
    pop <- 100 * popdisp * exp(growth) / (100 + popdisp * (exp(growth) - 1))
    mass[[yr]] <- data.frame(year = yr, mass_start = m_start, dispersing = P * m_start,
                             lost_kernel = P * m_start * (1 - w_sum), lost_domain = P * m_start * w_sum - sum(inside),
                             lost_unsuitable = lost_unsuit, mass_after_dispersal = sum(popdisp),
                             mass_after_growth = sum(pop))
    for (j in seq_along(levels)) yrs[[j]][is.infinite(yrs[[j]]) & pop > levels[j]] <- yr
    if (keep_density && yr %in% t) keep[[as.character(yr)]] <- pop
  }
  secs <- proc.time()[["elapsed"]] - t0
  tovec <- function(m) { v <- as.vector(t(m)); v[is.na(g$ei)] <- NA; v }
  list(density = if (keep_density) sapply(as.character(t), function(kk) tovec(keep[[kk]])) else NULL,
       years = sapply(yrs, tovec), mass = do.call(rbind, mass), status = "ok",
       params = c(base_params, list(cutoff = cutoff, R_c_requested_km = R_c, R_c_km = R_used, k_eff = R_used / u,
                                    kernel = kernel, kernel_half_cells = c(n_r, n_x), kernel_weight_sum = w_sum,
                                    tail_mass_beyond_cutoff = if (is.finite(R_c)) kernel_tail_2dt(R_c, u, p) else 0,
                                    cutoff_capped = is.finite(R_c) && R_c > R_used, grid_diag_km = h * sqrt(nr^2 + nc^2),
                                    fft_dim = c(Nr, Nc), seconds = secs)))
}

#' Dispersal kernel (DK) model, FFT convolution on an equal-area grid
#'
#' SpatRaster wrapper of [core_dispk()]; see there for the method, cutoff and mass diagnostics.
#'
#' @param pra An `sg_pra` object on a projected, equal-area grid with square cells (e.g. EPSG:3035).
#' @param init Initial density, SpatRaster from [init_density()] (percent of K).
#' @inheritParams core_dispk
#' @return An object of class `sg_dispk`: `density` (SpatRaster per `t`, percent of K), `years_to` (SpatRaster, one
#'   layer per level: first year above it, `Inf` = not within the horizon), `summary` (per `t` x level), `mass`,
#'   `params`, `status`. Skipped runs (`u` below `min_u_cells`) give a warning and `status = "spread below grid"`.
#' @export
spread_dispk <- function(pra, init, lmax, u, p, t = c(5, 10, 20, 50), P = 1, levels = c(0.1, 1, 10),
                         cutoff = c("auto", "k", "none"), k = 10, tail_safety = 0.1, kernel = c("point", "cell"),
                         cell_sub = 5, min_u_cells = 2) {
  if (terra::is.lonlat(pra$rast)) stop("spread_dispk() needs a projected equal-area grid; use spread_dispk_reference() on lon/lat")
  g <- as_grid(pra)
  res <- core_dispk(g, terra::values(init, mat = FALSE), lmax, u, p, t, P, levels, match.arg(cutoff), k, tail_safety,
                    match.arg(kernel), cell_sub, min_u_cells)
  if (!identical(res$status, "ok")) {
    warning(sprintf("u = %g km is %.2f cells (< %g): DK skipped (spread below grid)", u, u / g$cell_km, min_u_cells))
    return(structure(list(status = res$status, params = res$params), class = "sg_dispk"))
  }
  dispk_result(res, g, t)
}

#' Time-to-event over a target area from a per-cell years raster
#'
#' @param years SpatRaster of years (e.g. a layer of `years_to` of [spread_dispk()], or `years_to_reach` of
#'   [spread_radial()]; `Inf` = never).
#' @param pra The `sg_pra` object the raster belongs to.
#' @param target Optional SpatRaster on the PRA grid, 1 (or TRUE) inside the target area. NULL = the whole PRA.
#' @param pct Shares (percent) of the target's suitable area.
#' @param horizon Years simulated; events after it are censored (`Inf`, "never within horizon").
#' @return Data frame: `first_year` (first suitable target cell), `year_pct_<x>` (first year at which at least x % of
#'   the target's suitable area is reached), `censored_*` flags.
#' @export
time_to_event <- function(years, pra, target = NULL, pct = c(10, 50), horizon = Inf) {
  g <- as_grid(pra, target)
  tte_vec(terra::values(years, mat = FALSE), g, pct, horizon)
}

#' @export
print.sg_dispk <- function(x, ...) {
  pr <- x$params
  if (!identical(x$status, "ok")) {
    cat(sprintf("<sg_dispk> skipped: %s (u = %g km, cell %g km)\n", x$status, pr$u, pr$cell_km))
    return(invisible(x))
  }
  cat(sprintf("<sg_dispk> %s, lmax = %g, u = %g km, p = %g, P = %g, horizon %d yr, levels %s %% of K\n", pr$impl,
              pr$lmax, pr$u, pr$p, pr$P, pr$horizon, paste(pr$levels, collapse = " / ")))
  if (pr$impl == "fft") {
    cat(sprintf("  cutoff %s: R_c = %.0f km (%.1f u), kernel weight sum %.4f, tail mass beyond R_c %.2e, %.1f s\n",
                pr$cutoff, pr$R_c_km, pr$k_eff, pr$kernel_weight_sum, pr$tail_mass_beyond_cutoff, pr$seconds))
  }
  print(x$summary, row.names = FALSE, digits = 4)
  invisible(x)
}

# ---- internal ----

check_kernel_params <- function(lmax, p, P) {
  stopifnot(lmax >= 1, p > 1, P >= 0, P <= 1)
  if (p < 2) warning("p < 2: below the tutorial range (fat, near-Cauchy tail)")
}

check_levels <- function(levels) {
  stopifnot(is.numeric(levels), length(levels) >= 1, all(levels > 0), all(levels < 100))
}

# Core result (vectors) -> sg_dispk with SpatRasters
dispk_result <- function(res, g, t) {
  lv <- res$params$levels
  structure(list(density = grid_to_rast(g, res$density, paste0("t", t)),
                 years_to = grid_to_rast(g, res$years, paste0("years_to_", lv)),
                 summary = level_summary(res$density, g, t, lv), mass = res$mass, params = res$params,
                 status = "ok"),
            class = "sg_dispk")
}

# Smallest r with 100 * f(r) * h^2 <= target (percent of K)
auto_cutoff <- function(u, p, h, target) {
  c0 <- exp(lgamma((p + 1) / 2) - lgamma((p - 1) / 2)) / (u^2 * p * pi)
  ratio <- 100 * c0 * h^2 / target
  if (ratio <= 1) return(3 * h)
  max(3 * h, u * sqrt(p * (ratio^(2 / (p + 1)) - 1)))
}

# Kernel density values (km-2) on row offsets -n_r..n_r and column offsets -n_x..n_x, zero beyond R_c
kernel_matrix <- function(u, p, h, n_r, n_x, R_c, kernel, cell_sub) {
  dx <- matrix((-n_x):n_x, 2 * n_r + 1, 2 * n_x + 1, byrow = TRUE)
  dy <- matrix((-n_r):n_r, 2 * n_r + 1, 2 * n_x + 1)
  rr <- h * sqrt(dx^2 + dy^2)
  w <- if (kernel == "point") {
    kernel_2dt(rr, u, p)
  } else {
    s <- ((seq_len(cell_sub) - 0.5) / cell_sub - 0.5) * h
    acc <- 0
    for (sx in s) for (sy in s) acc <- acc + kernel_2dt(sqrt((h * dx + sx)^2 + (h * dy + sy)^2), u, p)
    acc / cell_sub^2
  }
  w[rr > R_c] <- 0
  w
}

# Place kernel offsets into an Nr x Nc matrix with the origin at [1, 1] (circular layout)
embed_kernel <- function(w, n_r, n_x, Nr, Nc) {
  K <- matrix(0, Nr, Nc)
  K[(((-n_r):n_r) %% Nr) + 1, (((-n_x):n_x) %% Nc) + 1] <- w
  K
}
