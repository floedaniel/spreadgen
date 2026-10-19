# Plain-data grid (no terra objects): the form in which model cores and Monte Carlo runs receive the PRA.
# Cell order is terra's (row-major from the top-left), so vectors convert to and from SpatRaster without reordering.

#' Convert a PRA object to a plain-data grid
#'
#' Model cores ([core_slg()], [core_radial_years()], [core_dispk()]) and [mc_spread()] work on this object only, so
#' nothing from terra crosses process boundaries in parallel runs. Convert back with [grid_to_rast()].
#'
#' @param pra An `sg_pra` object.
#' @param target Optional SpatRaster on the PRA grid, 1 / TRUE inside the target area (reporting region).
#' @return An object of class `sg_grid`: `nrow`, `ncol`, `ext` (xmin, xmax, ymin, ymax), `res` (x, y in CRS units),
#'   `cell_km` (cell side in km on projected grids, NA on lon/lat), `crs` (WKT), `lonlat`, and vectors over all cells:
#'   `ei` (0/1, NA = no data), `gi` (0-100), `area_km2`, `target` (logical, all TRUE when `target` is NULL).
#' @export
as_grid <- function(pra, target = NULL) {
  stopifnot(inherits(pra, "sg_pra"))
  r <- pra$rast
  v <- terra::values(r, mat = TRUE)
  e <- as.vector(terra::ext(r))
  tgt <- if (is.null(target)) rep(TRUE, terra::ncell(r)) else {
    g <- terra::values(target, mat = FALSE)
    !is.na(g) & g > 0
  }
  lonlat <- terra::is.lonlat(r)
  res <- terra::res(r)
  structure(list(nrow = terra::nrow(r), ncol = terra::ncol(r), ext = unname(e), res = res,
                 cell_km = if (lonlat) NA_real_ else res[1] / 1000, crs = terra::crs(r), lonlat = lonlat,
                 ei = unname(v[, "ei"]), gi = unname(v[, "gi"]), area_km2 = unname(v[, "area_km2"]), target = tgt,
                 settings = pra$settings),
            class = "sg_grid")
}

#' Values on a plain-data grid to SpatRaster (I/O edge)
#'
#' @param g An `sg_grid` object.
#' @param values Vector (one layer) or matrix (one column per layer) in cell order.
#' @param names Layer names.
#' @return SpatRaster, masked to cells with data.
#' @export
grid_to_rast <- function(g, values, names = NULL) {
  r <- terra::rast(nrows = g$nrow, ncols = g$ncol, xmin = g$ext[1], xmax = g$ext[2], ymin = g$ext[3],
                   ymax = g$ext[4], crs = g$crs, nlyrs = NCOL(values))
  values <- as.matrix(values)
  values[is.na(g$ei), ] <- NA
  terra::values(r) <- values
  if (!is.null(names)) names(r) <- names
  r
}

#' Cell-centre coordinates of a plain-data grid
#'
#' @param g An `sg_grid` object.
#' @param cells Cell numbers (default: all).
#' @return Two-column matrix (x, y) in the grid CRS (degrees or metres).
#' @export
grid_xy <- function(g, cells = seq_len(g$nrow * g$ncol)) {
  row <- (cells - 1) %/% g$ncol + 1
  col <- (cells - 1) %% g$ncol + 1
  cbind(x = g$ext[1] + (col - 0.5) * g$res[1], y = g$ext[4] - (row - 0.5) * g$res[2])
}

#' Distance (km) from every cell to the nearest of a set of cells, plain arithmetic
#'
#' @param g An `sg_grid` object.
#' @param from Cell numbers of the origins (e.g. entry cells).
#' @param method `"auto"` (lon/lat: `"sphere"`, projected: `"euclid"`), `"sphere"`, `"v21"` or `"euclid"`.
#' @return Numeric vector over all cells (NA where the grid has no data). Cost: length(from) x cells; for thousands
#'   of origins compute the distance once at the I/O edge (e.g. [spread_radial()] with `distance = "grid"`) and pass
#'   it to [mc_spread()].
#' @export
grid_dist_km <- function(g, from, method = c("auto", "sphere", "v21", "euclid")) {
  method <- match.arg(method)
  if (method == "auto") method <- if (g$lonlat) "sphere" else "euclid"
  if (method != "euclid" && !g$lonlat) stop("method '", method, "' needs a lon/lat grid")
  ok <- which(!is.na(g$ei))
  xy <- grid_xy(g)
  if (!g$lonlat) xy <- xy / 1000
  d <- rep(NA_real_, nrow(xy))
  digits <- if (method == "v21") 7 else NULL
  d[ok] <- Reduce(pmin, lapply(unique(from), function(i) {
    if (method == "euclid") sqrt((xy[ok, 1] - xy[i, 1])^2 + (xy[ok, 2] - xy[i, 2])^2)
    else dist_v21(xy[i, 1], xy[i, 2], xy[ok, 1], xy[ok, 2], digits)
  }))
  d
}

#' @export
print.sg_grid <- function(x, ...) {
  suit <- !is.na(x$ei) & x$ei == 1
  cat(sprintf("<sg_grid> %d x %d (%s), %d cells with data, %d suitable, %d suitable in target\n", x$nrow, x$ncol,
              if (x$lonlat) "lon/lat" else sprintf("%g km", x$cell_km), sum(!is.na(x$ei)), sum(suit),
              sum(suit & x$target)))
  invisible(x)
}

# ---- internal: summaries on plain vectors ----

# Long summary of a density matrix (cells x t, percent of K) at density levels (percent of K, strict >)
level_summary <- function(dens, g, t, levels) {
  dens <- as.matrix(dens)
  ok <- !is.na(g$ei)
  suit <- ok & g$ei == 1
  a <- g$area_km2
  a_suit <- sum(a[suit])
  do.call(rbind, lapply(seq_along(t), function(i) {
    d <- dens[, i]
    do.call(rbind, lapply(levels, function(lv) {
      hit <- ok & !is.na(d) & d > lv
      data.frame(t = t[i], level = lv, n_suitable = sum(suit), area_suitable_km2 = a_suit, n_cells = sum(hit),
                 area_km2 = sum(a[hit]), pct_suitable = 100 * sum(a[hit]) / a_suit,
                 mean_density_suitable = mean(d[suit]))
    }))
  }))
}

# Time-to-event on vectors: years per cell (Inf = never), over suitable target cells
tte_vec <- function(y, g, pct = c(10, 50), horizon = Inf) {
  sel <- g$target & !is.na(g$ei) & g$ei == 1
  y <- y[sel]
  a <- g$area_km2[sel]
  y[!is.na(y) & y > horizon] <- Inf
  out <- data.frame(n_suitable_target = sum(sel), area_suitable_target_km2 = sum(a),
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
