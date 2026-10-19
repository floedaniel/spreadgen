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
#' @param reach,abundant Density thresholds (percent of K) for `years_to_reach` / `years_to_abundant`.
#' @return An object of class `sg_dispk` (see [spread_dispk()]).
#' @export
spread_dispk_reference <- function(pra, init, lmax, u, p, t = 1:10, P = 1, distance = c("auto", "v21", "sphere", "euclid"),
                                   sphere_correction = NULL, reach = 1e-6, abundant = 1) {
  distance <- match.arg(distance)
  r <- pra$rast
  lonlat <- terra::is.lonlat(r)
  if (distance == "auto") distance <- if (lonlat) "sphere" else "euclid"
  if (distance != "euclid" && !lonlat) stop("distance '", distance, "' needs a lon/lat grid")
  if (is.null(sphere_correction)) sphere_correction <- lonlat
  check_kernel_params(lmax, p, P)
  v <- terra::values(r, mat = TRUE)
  ok <- which(!is.na(v[, "ei"]))
  xy <- terra::xyFromCell(r, ok)
  if (!lonlat) xy <- xy / 1000
  area <- v[ok, "area_km2"]
  suit <- v[ok, "ei"] == 1
  growth <- v[ok, "gi"] * log(lmax) / 100
  pop <- terra::values(init, mat = FALSE)[ok]
  pop[is.na(pop)] <- 0
  sources <- which(suit | pop > 0)
  digits <- if (distance == "v21") 7 else NULL
  t_max <- max(t)
  keep <- list()
  yr_reach <- yr_abund <- rep(Inf, length(ok))
  yr_reach[pop > reach] <- 0
  yr_abund[pop > abundant] <- 0
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
    yr_reach[is.infinite(yr_reach) & pop > reach] <- k
    yr_abund[is.infinite(yr_abund) & pop > abundant] <- k
    if (k %in% t) keep[[as.character(k)]] <- pop
  }
  dens <- terra::rast(lapply(as.character(t), function(k) cells_to_rast(r, ok, keep[[k]])))
  names(dens) <- paste0("t", t)
  new_dispk(dens, cells_to_rast(r, ok, yr_reach), cells_to_rast(r, ok, yr_abund), r, t, reach, abundant,
            mass = do.call(rbind, mass),
            params = list(model = "DK", impl = "reference", lmax = lmax, u = u, p = p, P = P, distance = distance,
                          sphere_correction = sphere_correction, horizon = t_max))
}

#' Dispersal kernel (DK) model, FFT convolution on an equal-area grid
#'
#' One dispersal step is a 2-D convolution of the density with the kernel value at the centre distance,
#' `f(r)`, multiplied by the landing-cell area (as the original; no renormalisation; on an equal-area grid this is
#' `f(r) h^2`), computed by zero-padded FFT (no wrap-around). Mass dispersing into unsuitable / no-data cells,
#' beyond the grid, or beyond the kernel cutoff is lost and reported per year.
#'
#' The kernel cutoff `R_c`: `"auto"` (default) is the distance beyond which all suitable cells together, each fully
#' occupied (100 % of K), add less than `tail_safety * reach` to a landing cell,
#' i.e. `n_suitable * 100 f(R_c) h^2 <= tail_safety * reach` (conservative; for fat tails this is the whole grid);
#' `"k"` uses `R_c = k * u`; `"none"` uses the whole grid (no truncation inside it).
#'
#' @inheritParams spread_dispk_reference
#' @param pra An `sg_pra` object on a projected, equal-area grid with square cells (e.g. EPSG:3035).
#' @param cutoff `"auto"`, `"k"` or `"none"`.
#' @param k Cutoff factor for `cutoff = "k"`.
#' @param tail_safety Factor for `cutoff = "auto"` (default 0.1).
#' @param kernel `"point"` (default, original discretisation) or `"cell"` (kernel averaged over `cell_sub` x
#'   `cell_sub` points per landing cell; deviation, better mass balance when `u` is near the cell size).
#' @param cell_sub Sub-points per side for `kernel = "cell"`.
#' @param min_u_cells Minimum `u` in cells; below it the model is skipped (returns `status = "spread below grid"`).
#' @return An object of class `sg_dispk`: `density` (SpatRaster per `t`, percent of K), `years_to_reach`,
#'   `years_to_abundant` (SpatRaster, first year above the threshold; `Inf` = not within the horizon), `summary`
#'   (per `t`), `mass` (per year, in percent-of-K x cells: `mass_start`, `dispersing`, `lost_kernel` (1 - kernel
#'   weight sum: truncation and discretisation), `lost_domain` (landed outside the grid or in no-data cells),
#'   `lost_unsuitable` (landed in unsuitable cells with data), `mass_after_dispersal`, `mass_after_growth`),
#'   `params` (incl. cutoff, kernel weight sum, analytic tail mass), `status` (`"ok"` or the skip reason).
#' @export
spread_dispk <- function(pra, init, lmax, u, p, t = c(5, 10, 20, 50), P = 1, reach = 1e-6, abundant = 1,
                         cutoff = c("auto", "k", "none"), k = 10, tail_safety = 0.1, kernel = c("point", "cell"),
                         cell_sub = 5, min_u_cells = 2) {
  cutoff <- match.arg(cutoff)
  kernel <- match.arg(kernel)
  r <- pra$rast
  if (terra::is.lonlat(r)) stop("spread_dispk() needs a projected equal-area grid; use spread_dispk_reference() on lon/lat")
  res <- terra::res(r)
  if (abs(res[1] - res[2]) > 1e-6 * res[1]) stop("spread_dispk() needs square cells")
  h <- res[1] / 1000
  check_kernel_params(lmax, p, P)
  if (reach < 1e-8) {
    warning("reach below 1e-8 % of K is within FFT round-off; set to 1e-8")
    reach <- 1e-8
  }
  t_max <- max(t)
  base_params <- list(model = "DK", impl = "fft", lmax = lmax, u = u, p = p, P = P, cell_km = h, horizon = t_max)
  if (u / h < min_u_cells) {
    warning(sprintf("u = %g km is %.2f cells (< %g): DK skipped (spread below grid)", u, u / h, min_u_cells))
    return(structure(list(status = "spread below grid", params = base_params), class = "sg_dispk"))
  }
  nr <- terra::nrow(r)
  nc <- terra::ncol(r)
  R_c <- switch(cutoff,
                none = Inf,
                k = k * u,
                auto = auto_cutoff(u, p, h, tail_safety * reach / max(1, sum(terra::values(r[["ei"]], mat = FALSE) == 1, na.rm = TRUE))))
  diag_km <- h * sqrt(nr^2 + nc^2)
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

  ei <- terra::as.matrix(r[["ei"]], wide = TRUE)
  has_data <- !is.na(ei)
  suit <- has_data & ei == 1
  area <- terra::as.matrix(r[["area_km2"]], wide = TRUE)
  area[is.na(area)] <- 0
  growth <- terra::as.matrix(r[["gi"]], wide = TRUE) * log(lmax) / 100
  growth[is.na(growth)] <- 0
  pop <- terra::as.matrix(init, wide = TRUE)
  pop[is.na(pop)] <- 0
  yr_reach <- ifelse(pop > reach, 0, Inf)
  yr_abund <- ifelse(pop > abundant, 0, Inf)
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
    yr_reach[is.infinite(yr_reach) & pop > reach] <- yr
    yr_abund[is.infinite(yr_abund) & pop > abundant] <- yr
    if (yr %in% t) keep[[as.character(yr)]] <- pop
  }
  secs <- proc.time()[["elapsed"]] - t0
  to_r <- function(m) terra::mask(terra::setValues(r[["ei"]], as.vector(t(m))), r[["ei"]])
  dens <- terra::rast(lapply(as.character(t), function(kk) to_r(keep[[kk]])))
  names(dens) <- paste0("t", t)
  new_dispk(dens, to_r(yr_reach), to_r(yr_abund), r, t, reach, abundant, mass = do.call(rbind, mass),
            params = c(base_params, list(cutoff = cutoff, R_c_requested_km = R_c, R_c_km = R_used, k_eff = R_used / u,
                                         kernel = kernel, kernel_half_cells = c(n_r, n_x), kernel_weight_sum = w_sum,
                                         tail_mass_beyond_cutoff = if (is.finite(R_c)) kernel_tail_2dt(R_c, u, p) else 0,
                                         cutoff_capped = is.finite(R_c) && R_c > R_used, grid_diag_km = diag_km,
                                         fft_dim = c(Nr, Nc), seconds = secs)))
}

#' Time-to-event over a region from a per-cell years raster
#'
#' @param years SpatRaster of years (e.g. `years_to_reach` of [spread_dispk()] or [spread_radial()]; `Inf` = never).
#' @param pra The `sg_pra` object the raster belongs to.
#' @param region Optional SpatRaster on the PRA grid, 1 (or TRUE) inside the reporting region (e.g. Norway).
#'   NULL = the whole PRA.
#' @param pct Shares (percent) of the region's suitable area.
#' @param horizon Years simulated; events after it are censored (`Inf`, "never within horizon").
#' @return Data frame: `first_year` (first suitable region cell reached), `year_pct_<x>` (first year at which at
#'   least x % of the region's suitable area is reached), `censored_*` flags.
#' @export
time_to_event <- function(years, pra, region = NULL, pct = c(10, 50), horizon = Inf) {
  y <- terra::values(years, mat = FALSE)
  ei <- terra::values(pra$rast[["ei"]], mat = FALSE)
  a <- terra::values(pra$rast[["area_km2"]], mat = FALSE)
  inreg <- if (is.null(region)) rep(TRUE, length(y)) else { g <- terra::values(region, mat = FALSE); !is.na(g) & g > 0 }
  sel <- inreg & !is.na(ei) & ei == 1
  y <- y[sel]
  a <- a[sel]
  y[y > horizon] <- Inf
  out <- data.frame(n_suitable_region = sum(sel), area_suitable_region_km2 = sum(a),
                    first_year = if (length(y)) min(y) else Inf)
  o <- order(y)
  cum <- cumsum(a[o]) / sum(a)
  for (x in pct) {
    j <- which(cum >= x / 100 - 1e-12)[1]
    out[[paste0("year_pct_", x)]] <- if (is.na(j)) Inf else y[o][j]
  }
  for (nm in grep("year", names(out), value = TRUE)) out[[paste0("censored_", nm)]] <- is.infinite(out[[nm]])
  out
}

#' @export
print.sg_dispk <- function(x, ...) {
  pr <- x$params
  if (!identical(x$status, "ok")) {
    cat(sprintf("<sg_dispk> skipped: %s (u = %g km, cell %g km)\n", x$status, pr$u, pr$cell_km))
    return(invisible(x))
  }
  cat(sprintf("<sg_dispk> %s, lmax = %g, u = %g km, p = %g, P = %g, horizon %d yr\n", pr$impl, pr$lmax, pr$u, pr$p,
              pr$P, pr$horizon))
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

new_dispk <- function(dens, yr_reach, yr_abund, r, t, reach, abundant, mass, params) {
  names(yr_reach) <- "years_to_reach"
  names(yr_abund) <- "years_to_abundant"
  structure(list(density = dens, years_to_reach = yr_reach, years_to_abundant = yr_abund,
                 summary = density_summary(dens, r, t, reach, abundant), mass = mass,
                 params = c(params, list(reach = reach, abundant = abundant)), status = "ok"),
            class = "sg_dispk")
}

cells_to_rast <- function(r, ok, x) {
  v <- rep(NA_real_, terra::ncell(r))
  v[ok] <- x
  terra::setValues(r[["ei"]], v)
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
