# Model RR (paper: model B): radial range expansion; Rand-RR: random cells inside the radial range (hybrid A + B).

#' Radial range expansion (RR) model
#'
#' A suitable cell is invaded at time `t` when its distance to the nearest entry cell is at most `RR * t`
#' (union of discs around the entry cells only; newly invaded cells do not start new discs). The year a cell is
#' reached is `ceil(d / RR)` (entry cells: 0). Unsuitable cells are never invaded.
#'
#' @param pra An `sg_pra` object.
#' @param entry Entry points: a two-column matrix / data frame of coordinates (snapped with [snap_entry()] to the
#'   nearest suitable cell), or an integer vector of PRA cell numbers (must be suitable).
#' @param RR Radial rate of spread (km/yr), > 0.
#' @param t Years (vector) at which the invaded set is returned.
#' @param distance Distance to the nearest entry cell. `"sphere"` (default): great circle on a sphere of radius
#'   20000/pi km (as the original), cosine clamped to \[-1, 1\], no rounding. `"v21"`: the original formula with the
#'   cosine rounded to 7 digits before `acos` (exact regression). `"grid"`: `terra::distance()` (WGS84 geodesic on
#'   lon/lat grids, Euclidean on projected grids; fastest with many entry cells).
#' @param reach `"disc"` (default, Robinet et al. 2012): straight-line distance. `"land"`: shortest path through
#'   cells with data (`terra::gridDist()`), so sea and no-data gaps must be travelled around (not in the original).
#' @param entry_crs CRS of `entry` coordinates when it differs from the PRA grid.
#' @return An object of class `sg_radial`: `invaded` (SpatRaster, one 0/1 layer per `t`), `years_to_reach`
#'   (SpatRaster, `Inf` for unsuitable or unreachable cells), `distance_km` (SpatRaster, distance to the nearest
#'   entry cell), `entry` (cells), `summary` (data frame per `t`), `params`.
#' @export
spread_radial <- function(pra, entry, RR, t = c(5, 10, 20, 50), distance = c("sphere", "v21", "grid"),
                          reach = c("disc", "land"), entry_crs = NULL) {
  distance <- match.arg(distance)
  reach <- match.arg(reach)
  stopifnot(inherits(pra, "sg_pra"), length(RR) == 1, RR > 0, all(t >= 0))
  if (reach == "land" && distance != "sphere") stop("reach = 'land' uses its own grid distance; leave distance at its default")
  r <- pra$rast
  ent <- entry_cells(pra, entry, entry_crs)
  d <- entry_distance_km(r, ent, distance, reach)
  suit <- r[["ei"]] == 1
  # years to reach: inclusive test d <= RR * t  <=>  t >= ceil(d / RR) for integer t
  yrs <- terra::ifel(suit & is.finite(d), ceiling(d / RR), Inf)
  yrs <- terra::mask(yrs, r[["ei"]])
  names(yrs) <- "years_to_reach"
  inv <- terra::rast(lapply(t, function(tt) terra::ifel(suit & d <= RR * tt, 1, 0)))
  inv <- terra::mask(inv, r[["ei"]])
  names(inv) <- paste0("t", t)
  out <- list(invaded = inv, years_to_reach = yrs, distance_km = d, entry = ent,
              summary = occupancy_summary(inv, r, t, d),
              params = list(model = "RR", RR = RR, distance = distance, reach = reach, n_entry = length(ent)))
  class(out) <- "sg_radial"
  out
}

#' Random radial range expansion (Rand-RR)
#'
#' The radial range of [spread_radial()] gives the cells that can be invaded; the number invaded follows the
#' logistic of invaded cells, `round(N_t / 100 * n_suitable)` with `N_t = logistic_pct(n0, r, t)`, chosen at random
#' inside the range (all range cells when fewer). Cells are drawn in one random order per run, so the invaded set is
#' nested in time (the original drew independently for each call). Uses the caller's RNG stream (seeded by
#' `mc_spread()`); call `set.seed()` first for a reproducible single run.
#'
#' @inheritParams spread_radial
#' @param n0 Initial share of suitable cells invaded, percent (e.g. `n0_from(1, n_suitable)`).
#' @param r Relative rate of increase of invaded cells (yr-1), e.g. from [r_from_invaded()].
#' @return An object of class `sg_radial` (model `"Rand-RR"`).
#' @export
spread_radialrand <- function(pra, entry, RR, n0, r, t = c(5, 10, 20, 50),
                              distance = c("sphere", "v21", "grid"), reach = c("disc", "land"), entry_crs = NULL) {
  rad <- spread_radial(pra, entry, RR, t, match.arg(distance), match.arg(reach), entry_crs)
  ras <- pra$rast
  ei <- terra::values(ras[["ei"]], mat = FALSE)
  suit_cells <- which(ei == 1)
  n_suit <- length(suit_cells)
  prio <- stats::runif(length(ei))
  dv <- terra::values(rad$distance_km, mat = FALSE)
  layers <- lapply(seq_along(t), function(i) {
    in_range <- suit_cells[dv[suit_cells] <= RR * t[i]]
    n_inv <- round(logistic_pct(n0, r, t[i]) * n_suit / 100)
    keep <- if (n_inv >= length(in_range)) in_range else in_range[order(prio[in_range])[seq_len(n_inv)]]
    v <- ifelse(is.na(ei), NA, 0)
    v[keep] <- 1
    terra::setValues(ras[["ei"]], v)
  })
  inv <- terra::rast(layers)
  names(inv) <- paste0("t", t)
  rad$invaded <- inv
  rad$summary <- occupancy_summary(inv, ras, t, rad$distance_km)
  rad$params <- c(rad$params, list(model = "Rand-RR", n0 = n0, r = r))
  rad
}

#' @export
print.sg_radial <- function(x, ...) {
  cat(sprintf("<sg_radial> model %s, RR = %g km/yr, %d entry cell(s), distance = %s, reach = %s\n",
              x$params$model, x$params$RR, x$params$n_entry, x$params$distance, x$params$reach))
  print(x$summary, row.names = FALSE, digits = 4)
  invisible(x)
}

# ---- internal ----

# Entry cells (suitable) from coordinates or cell numbers
entry_cells <- function(pra, entry, entry_crs = NULL) {
  if (is.null(dim(entry))) {
    ent <- unique(as.integer(entry))
    ei <- terra::values(pra$rast[["ei"]], mat = FALSE)[ent]
    if (anyNA(ei) || any(ei != 1)) stop("entry cells must be suitable cells of the PRA")
    return(ent)
  }
  unique(snap_entry(pra, entry, crs = entry_crs)$cell)
}

# Distance (km) from every cell to the nearest entry cell
entry_distance_km <- function(r, ent, distance, reach) {
  tmpl <- r[["ei"]]
  if (reach == "land") {
    # cells with data are passable, NA (sea / no data) is a barrier; entry cells are the targets (value 0)
    g <- terra::ifel(is.na(tmpl), NA, 1)
    g[ent] <- 0
    d <- terra::gridDist(g, target = 0) / 1000
  } else if (distance == "grid") {
    e <- terra::rast(tmpl)
    e[ent] <- 1
    d <- terra::distance(e) / 1000
  } else {
    xy <- terra::xyFromCell(tmpl, seq_len(terra::ncell(tmpl)))
    if (!terra::is.lonlat(tmpl)) {
      xy <- terra::geom(terra::project(terra::vect(xy, crs = terra::crs(tmpl)), "EPSG:4326"))[, c("x", "y")]
    }
    ok <- which(!is.na(terra::values(tmpl, mat = FALSE)))
    dv <- rep(NA_real_, nrow(xy))
    digits <- if (distance == "v21") 7 else NULL
    dv[ok] <- Reduce(pmin, lapply(ent, function(i) dist_v21(xy[i, 1], xy[i, 2], xy[ok, 1], xy[ok, 2], digits)))
    d <- terra::setValues(tmpl, dv)
  }
  d <- terra::mask(d, tmpl)
  names(d) <- "distance_km"
  d
}

#' Great-circle distance as in Prog_SpreadModule_V21-web.R
#'
#' `R = (20000/pi) * acos(round(sin(y1) sin(y2) + cos(y1) cos(y2) cos(x2 - x1), 7))`, coordinates in degrees.
#' @param x1,y1 Longitude / latitude of the origin (degrees).
#' @param x2,y2 Longitude / latitude vectors (degrees).
#' @param round_digits Digits of the cosine rounding (7 in the original); `NULL` = no rounding, cosine clamped to
#'   \[-1, 1\].
#' @return Distance in km.
#' @export
dist_v21 <- function(x1, y1, x2, y2, round_digits = 7) {
  k <- pi / 180
  cs <- sin(k * y1) * sin(k * y2) + cos(k * y1) * cos(k * y2) * cos(k * (x2 - x1))
  cs <- if (is.null(round_digits)) pmin(1, pmax(-1, cs)) else round(cs, round_digits)
  20000 * acos(cs) / pi
}

# Per-layer summary of an occupancy stack (0/1)
occupancy_summary <- function(inv, r, t, d = NULL) {
  v <- terra::values(inv, mat = TRUE)
  ei <- terra::values(r[["ei"]], mat = FALSE)
  a <- terra::values(r[["area_km2"]], mat = FALSE)
  dv <- if (is.null(d)) NULL else terra::values(d, mat = FALSE)
  suit <- !is.na(ei) & ei == 1
  n_suit <- sum(suit)
  a_suit <- sum(a[suit])
  do.call(rbind, lapply(seq_along(t), function(i) {
    hit <- !is.na(v[, i]) & v[, i] == 1
    data.frame(t = t[i], n_suitable = n_suit, area_suitable_km2 = a_suit, n_invaded = sum(hit),
               area_invaded_km2 = sum(a[hit]), pct_suitable_invaded = 100 * sum(a[hit]) / a_suit,
               pniche = 100 * sum(hit) / n_suit,
               max_distance_km = if (!is.null(dv) && any(hit)) max(dv[hit]) else NA_real_)
  }))
}
