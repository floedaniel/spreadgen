# Shared test inputs: the V21 D. virgifera baseline and small synthetic grids.

v21 <- function() readRDS(test_path("fixtures", "v21_dvv_baseline.rds"))

# Original D. virgifera inputs (downloaded by data-raw/download_original.R; tests skip when absent)
dvv_file <- function(name) {
  f <- file.path(test_path(), "..", "..", "inst", "original", "Code for case studies",
                 "1-Diabrotica virgifera virgifera", name)
  if (!file.exists(f)) f <- system.file("original", "Code for case studies", "1-Diabrotica virgifera virgifera",
                                        name, package = "spreadgen")
  if (!nzchar(f) || !file.exists(f)) testthat::skip(paste("original input not downloaded:", name))
  f
}

# PRA from the fixture's CLIMEX table (EI / GI already habitat-masked and rescaled)
pra_v21 <- function(b = v21()) prep_climex(b$tab, rescale_gi = FALSE, area_method = "v21")

# Projected 1 km test grid (EPSG:3035), suitability increasing west -> east
synth_suit <- function(nr = 20, nc = 30, res_m = 1000) {
  r <- terra::rast(nrows = nr, ncols = nc, xmin = 4e6, xmax = 4e6 + nc * res_m, ymin = 3e6, ymax = 3e6 + nr * res_m,
                   crs = "EPSG:3035")
  terra::values(r) <- rep(seq(0, 1, length.out = nc), times = nr)
  names(r) <- "suit"
  r
}
