test_that("carrying capacity and N0 match the D. virgifera example", {
  K <- k_carrying(1579, 0.2, 200e6)
  expect_equal(K, 6.316e10)
  expect_equal(signif(K, 2), 6.3e10)
  expect_equal(n0_from(100, 6.3e10), 1.587e-7, tolerance = 1e-3)
  expect_equal(signif(n0_from(100, 6.3e10), 2), 1.6e-7)
})

test_that("r_from_invaded reproduces r = 0.33 (D. virgifera, model A)", {
  expect_equal(round(r_from_invaded(358, 3104, 1, 18), 2), 0.33)
  expect_error(r_from_invaded(3104, 3104, 1, 18))
})

test_that("prep_suitability: threshold, GI mapping and rescaling", {
  s <- synth_suit()
  p <- prep_suitability(s, threshold = 0.5)
  v <- terra::values(p$rast, mat = TRUE)
  sv <- terra::values(s, mat = FALSE)
  expect_s3_class(p, "sg_pra")
  expect_equal(v[, "ei"], as.numeric(sv >= 0.5))
  expect_equal(v[, "gi"], 100 * sv * (sv >= 0.5))
  expect_equal(v[, "area_km2"], rep(1, length(sv)))

  s2 <- s * 0.5
  expect_equal(max(terra::values(prep_suitability(s2, 0.1)$rast[["gi"]])), 50)
  expect_equal(max(terra::values(prep_suitability(s2, 0.1, gi_rescale = "domain")$rast[["gi"]])), 100)
  expect_equal(max(terra::values(prep_suitability(s, 0.1, gi_fun = function(x) 50 * x)$rast[["gi"]])), 50)
  expect_error(prep_suitability(s * 2))
})

test_that("prep_suitability: habitat, elevation and NA masks", {
  s <- synth_suit()
  s[1:5] <- NA
  h <- s * 0 + 1
  h[1, ] <- 0                                   # first row: no host
  el <- s * 0 + 100
  el[2, ] <- 2000                               # second row: too high
  p <- prep_suitability(s, 0, habitat = h, elev = el, elevmax = 1000)
  ei <- terra::values(p$rast[["ei"]], mat = FALSE)
  expect_true(all(is.na(ei[1:5])))
  expect_true(all(ei[terra::cellFromRowCol(s, 1, 6:30)] == 0))
  expect_true(all(ei[terra::cellFromRowCol(s, 2, 1:30)] == 0))
  expect_true(all(ei[terra::cellFromRowCol(s, 3, 1:30)] == 1))
  expect_true(all(terra::values(p$rast[["gi"]], mat = FALSE)[ei == 0] == 0, na.rm = TRUE))
})

test_that("prep_climex on the original inputs reproduces the V21 point table", {
  b <- v21()
  p <- prep_climex(dvv_file("ClimexOutput.txt"), habitat = terra::rast(dvv_file("habitat.tif")))
  cells <- terra::cellFromXY(p$rast, cbind(b$tab$lon, b$tab$lat))
  v <- terra::values(p$rast, mat = TRUE)[cells, ]
  expect_equal(v[, "ei"], as.numeric(b$tab$EI > 0))
  suit <- b$tab$EI > 0                          # spreadgen sets GI = 0 outside suitable cells (deviation D3)
  expect_equal(v[suit, "gi"], b$tab$GI[suit], tolerance = 1e-12)
  expect_true(all(v[!suit, "gi"] == 0))
  expect_equal(sum(v[, "ei"]), 3104)
  expect_equal(sum(!is.na(terra::values(p$rast[["ei"]]))), nrow(b$tab))
})

test_that("cell_area: v21 method reproduces the original areas; terra (ellipsoid) within 1.5 %", {
  b <- v21()
  p <- pra_v21(b)
  g <- b$grid
  cells <- terra::cellFromXY(p$rast, cbind(b$tab$lon, b$tab$lat))
  key <- function(lat, lon) paste(round(lat, 3), round(lon, 3))
  orig <- g$area_km2[match(key(b$tab$lat, b$tab$lon), key(g$lat, g$lon))]
  a21 <- terra::values(p$rast[["area_km2"]], mat = FALSE)[cells]
  expect_equal(a21, orig, tolerance = 1e-3)    # original rounds dX to 0.01 degree
  a_terra <- terra::values(cell_area(p$rast[["ei"]], "terra"), mat = FALSE)[cells]
  expect_lt(max(abs(a_terra / a21 - 1)), 0.015)   # sphere (40,000 km circumference) vs WGS84 ellipsoid
  expect_error(cell_area(synth_suit(), "v21"))
})

test_that("snap_entry moves points to the nearest suitable cell", {
  s <- synth_suit()
  p <- prep_suitability(s, threshold = 0.5)       # suitable from column 15 eastwards (x >= 4,014,000 m)
  sn <- snap_entry(p, cbind(4002500, 3010500))    # column 3, unsuitable
  expect_equal(terra::colFromCell(p$rast, sn$cell), 16)
  expect_equal(sn$dist_km, 13, tolerance = 1e-9)
  sn2 <- snap_entry(p, cbind(4020500, 3010500))   # already suitable: stays
  expect_equal(sn2$dist_km, 0)
  expect_equal(snap_entry(p, cbind(4002500, 3010500), to = "any")$dist_km, 0)
  # lon/lat input reprojected to the grid CRS
  ll <- terra::geom(terra::project(terra::vect(cbind(4020500, 3010500), crs = "EPSG:3035"), "EPSG:4326"))[, c("x", "y")]
  expect_equal(snap_entry(p, matrix(ll, 1), crs = "EPSG:4326")$cell, sn2$cell)
})

test_that("snap_entry on the V21 grid: Belgrade is a suitable cell", {
  p <- pra_v21()
  sn <- snap_entry(p, cbind(20.30, 44.82))
  expect_equal(unname(unlist(sn[, c("x_cell", "y_cell")])), c(20.25, 44.75))
  expect_lt(sn$dist_km, 15)
})
