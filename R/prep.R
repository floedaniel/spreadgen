# Preparation of the PRA (pest risk analysis area) object: suitable mask (EI proxy), growth index (GI proxy), cell area.

#' Prepare a PRA object from a suitability raster
#'
#' Converts an SDM / climate-match suitability raster into the inputs of the spread models: a suitable-cell mask
#' (`ei`, the proxy for CLIMEX EI > 0), a growth index 0-100 (`gi`, the proxy for the rescaled CLIMEX GI) and the
#' cell area. Mapping suitability to GI is an assumption (CLIMEX GI is a population growth index; suitability is
#' not); change it with `gi_fun`.
#'
#' @param suit SpatRaster, one layer, suitability in \[0, 1\]. NA = no data (never suitable).
#' @param threshold Numeric, suitability at or above which a cell is suitable (default 0.1).
#' @param gi_rescale `"none"` (default): GI = `gi_fun(suit)`. `"domain"`: GI rescaled so that its maximum over all
#'   cells with data is 100, as the original code did for CLIMEX GI (results then depend on the extent).
#' @param gi_fun Function mapping suitability to GI (0-100). Default `function(s) 100 * s`.
#' @param habitat Optional SpatRaster (same geometry), host presence or proportion; cells with habitat <= 0 or NA
#'   are unsuitable.
#' @param elev Optional SpatRaster (same geometry), elevation in m.
#' @param elevmax Numeric, maximum elevation (m) for establishment; cells at or above it (or with NA elevation)
#'   are unsuitable. Ignored when `elev` is NULL.
#' @param area_method Passed to [cell_area()].
#' @return An object of class `sg_pra`: a list with `rast` (SpatRaster with layers `ei` (0/1), `gi` (0-100, 0 where
#'   unsuitable) and `area_km2`), and the settings used.
#' @export
prep_suitability <- function(suit, threshold = 0.1, gi_rescale = c("none", "domain"), gi_fun = function(s) 100 * s,
                             habitat = NULL, elev = NULL, elevmax = NULL, area_method = c("terra", "v21")) {
  gi_rescale <- match.arg(gi_rescale)
  stopifnot(inherits(suit, "SpatRaster"), terra::nlyr(suit) == 1, is.numeric(threshold), length(threshold) == 1)
  rng <- terra::global(suit, range, na.rm = TRUE)
  if (rng[1, 1] < 0 || rng[1, 2] > 1) stop("suit must be in [0, 1]")
  gi <- gi_fun(suit)
  if (gi_rescale == "domain") gi <- gi * 100 / terra::global(gi, "max", na.rm = TRUE)[1, 1]
  ei <- suit >= threshold
  new_pra(ei, gi, habitat, elev, elevmax, area_method = match.arg(area_method),
          settings = list(source = "suitability", threshold = threshold, gi_rescale = gi_rescale))
}

#' Prepare a PRA object from CLIMEX output (legacy / regression path)
#'
#' Reads the point table used by the original code (latitude, longitude, EI, GI; tab-delimited, no header) onto a
#' regular lon/lat grid, rescales GI to a maximum of 100 over the table and applies the optional habitat and
#' elevation masks as `Prog_SpreadModule_V21-web.R` does.
#'
#' @param x Path to `ClimexOutput.txt`, or a data frame with columns `lat`, `lon`, `EI`, `GI`.
#' @param res Grid resolution in degrees (default 0.5).
#' @param habitat Optional SpatRaster (lon/lat); values extracted at the points, NA treated as 0 (as the original).
#' @param elev,elevmax Optional elevation SpatRaster (lon/lat) and maximum elevation (m).
#' @param area_method Passed to [cell_area()]; `"v21"` reproduces the original cell areas.
#' @param rescale_gi Rescale GI to max 100 over the table (TRUE, as the original). FALSE when GI is already 0-100.
#' @return An `sg_pra` object on a lon/lat grid (EPSG:4326). Cells without a CLIMEX point are NA.
#' @export
prep_climex <- function(x, res = 0.5, habitat = NULL, elev = NULL, elevmax = NULL, area_method = c("v21", "terra"),
                        rescale_gi = TRUE) {
  tab <- if (is.character(x)) {
    stats::setNames(utils::read.table(x, header = FALSE)[, 1:4], c("lat", "lon", "EI", "GI"))
  } else as.data.frame(x)[, c("lat", "lon", "EI", "GI")]
  tab$EI[is.na(tab$EI)] <- 0
  if (rescale_gi) tab$GI <- tab$GI * 100 / max(tab$GI, na.rm = TRUE)
  tab$GI[is.na(tab$GI)] <- 0
  xy <- cbind(tab$lon, tab$lat)
  if (!is.null(habitat)) {
    h <- terra::extract(habitat, xy)[, 1]
    h[is.na(h)] <- 0
    tab$EI <- tab$EI * (h > 0)
    tab$GI <- tab$GI * (h > 0)
  }
  if (!is.null(elev) && !is.null(elevmax)) {
    e <- terra::extract(elev, xy)[, 1]
    ok <- !is.na(e) & e < elevmax
    tab$EI <- tab$EI * ok
    tab$GI <- tab$GI * ok
  }
  tmpl <- grid_from_points(tab$lon, tab$lat, res)
  cells <- terra::cellFromXY(tmpl, xy)
  ei <- gi <- tmpl
  ei[cells] <- as.numeric(tab$EI > 0)
  gi[cells] <- tab$GI
  new_pra(ei, gi, area_method = match.arg(area_method),
          settings = list(source = "climex", threshold = NA_real_, gi_rescale = "domain"))
}

#' Cell area in km2
#'
#' @param x SpatRaster.
#' @param method `"terra"` (default): `terra::cellSize()` (exact on the ellipsoid for lon/lat, planar for projected
#'   grids). `"v21"`: the original approximation for lon/lat grids,
#'   `(40000/360)^2 * dx * dy * cos(lat)` km2 (sphere with a 40,000 km circumference).
#' @return SpatRaster `area_km2`, NA where `x` is NA.
#' @export
cell_area <- function(x, method = c("terra", "v21")) {
  method <- match.arg(method)
  if (method == "terra") {
    a <- terra::cellSize(x[[1]], unit = "km", mask = TRUE)
  } else {
    if (!terra::is.lonlat(x)) stop("method = 'v21' needs a lon/lat grid")
    r <- terra::res(x)
    a <- terra::init(x[[1]], "y")
    a <- (40000 / 360)^2 * r[1] * r[2] * abs(cos(a * pi / 180))
    a <- terra::mask(a, x[[1]])
  }
  names(a) <- "area_km2"
  a
}

#' Carrying capacity per cell
#'
#' @param area_km2 Cell area (km2); numeric or SpatRaster.
#' @param host_prop Proportion of the cell covered by host (0-1); numeric or SpatRaster (host-weighted K).
#' @param density_km2 Maximum density per km2 of host.
#' @return K, individuals (or units) per cell. Example (D. virgifera): `k_carrying(1579, 0.2, 200e6)` = 6.3e10.
#' @export
k_carrying <- function(area_km2, host_prop, density_km2) area_km2 * host_prop * density_km2

#' Initial density in percent of K
#'
#' @param n0 Number of founders (or invaded cells) at t = 0.
#' @param K Carrying capacity (or number of suitable cells).
#' @return `100 * n0 / K` (percent). Example: `n0_from(100, 6.3e10)` = 1.6e-7 %.
#' @export
n0_from <- function(n0, K) 100 * n0 / K

#' Snap entry points to the nearest suitable cell
#'
#' Points outside the suitable area move to the nearest suitable cell (great-circle distance on lon/lat grids,
#' Euclidean on projected grids; the original used squared degree differences). Ties: first cell in raster order.
#'
#' @param pra An `sg_pra` object.
#' @param xy Two-column matrix or data frame of coordinates (x/lon, y/lat).
#' @param crs CRS of `xy`; default: the CRS of the PRA grid.
#' @param to `"suitable"` (default) or `"any"` (any cell with data, as the original did for presence files).
#' @return Data frame: input coordinates, `cell`, `x_cell`, `y_cell` (cell centre, PRA CRS) and `dist_km`.
#' @export
snap_entry <- function(pra, xy, crs = NULL, to = c("suitable", "any")) {
  to <- match.arg(to)
  r <- pra$rast
  xy <- as.matrix(xy)[, 1:2, drop = FALSE]
  pts <- xy
  if (!is.null(crs)) {
    pts <- terra::geom(terra::project(terra::vect(xy, crs = crs), terra::crs(r)))[, c("x", "y"), drop = FALSE]
  }
  target <- if (to == "suitable") which(terra::values(r[["ei"]], mat = FALSE) == 1) else which(!is.na(terra::values(r[["ei"]], mat = FALSE)))
  if (!length(target)) stop("no target cells")
  txy <- terra::xyFromCell(r, target)
  lonlat <- terra::is.lonlat(r)
  out <- lapply(seq_len(nrow(pts)), function(i) {
    d <- terra::distance(pts[i, , drop = FALSE], txy, lonlat = lonlat)[1, ]
    j <- which.min(d)
    data.frame(x = xy[i, 1], y = xy[i, 2], cell = target[j], x_cell = txy[j, 1], y_cell = txy[j, 2], dist_km = d[j] / 1000)
  })
  do.call(rbind, out)
}

#' @export
print.sg_pra <- function(x, ...) {
  v <- terra::values(x$rast, mat = TRUE)
  ok <- !is.na(v[, "ei"])
  cat(sprintf("<sg_pra> %d x %d cells (%s), res %s, source: %s\n", terra::nrow(x$rast), terra::ncol(x$rast),
              if (terra::is.lonlat(x$rast)) "lon/lat" else "projected", paste(signif(terra::res(x$rast), 4), collapse = " x "),
              x$settings$source))
  cat(sprintf("  cells with data %d, suitable %d (%.0f km2), GI range %.1f-%.1f, gi_rescale = %s\n",
              sum(ok), sum(v[ok, "ei"] == 1), sum(v[ok & v[, "ei"] == 1, "area_km2"]),
              min(v[ok, "gi"]), max(v[ok, "gi"]), x$settings$gi_rescale))
  invisible(x)
}

# ---- internal ----

new_pra <- function(ei, gi, habitat = NULL, elev = NULL, elevmax = NULL, area_method = "terra", settings = list()) {
  na <- is.na(ei) | is.na(gi)
  if (!is.null(habitat)) {
    h <- terra::subst(habitat, NA, 0)
    ei <- ei * (h > 0)
  }
  if (!is.null(elev) && !is.null(elevmax)) {
    ok <- terra::subst(elev < elevmax, NA, 0)
    ei <- ei * ok
  }
  ei <- terra::mask(ei, na, maskvalues = 1)
  gi <- terra::mask(gi * (ei > 0), na, maskvalues = 1)
  r <- c(ei, gi, cell_area(ei, area_method))
  names(r) <- c("ei", "gi", "area_km2")
  structure(list(rast = r, settings = c(settings, list(area_method = area_method,
                                                       habitat = !is.null(habitat), elevmax = elevmax))),
            class = "sg_pra")
}

# Regular lon/lat grid whose cell centres are the given points
grid_from_points <- function(lon, lat, res) {
  e <- terra::ext(min(lon) - res / 2, max(lon) + res / 2, min(lat) - res / 2, max(lat) + res / 2)
  terra::rast(e, resolution = res, crs = "EPSG:4326")
}

# Cell values of a PRA as a data frame (cells with data only)
pra_cells <- function(pra) {
  v <- terra::values(pra$rast, mat = TRUE)
  cell <- which(!is.na(v[, "ei"]))
  data.frame(cell = cell, ei = v[cell, "ei"], gi = v[cell, "gi"], area_km2 = v[cell, "area_km2"])
}
