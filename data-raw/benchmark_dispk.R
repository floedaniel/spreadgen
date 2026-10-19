# M3 benchmark of the FFT dispersal-kernel model (spread_dispk) on the production-type grid:
# land in Europe north of 45 N, EPSG:3035, at 5 / 10 / 20 km (same construction as the sibling scripts' land grid).
# Runtime does not depend on the species (FFT size = grid + kernel), so suitability is synthetic (0.5 on land).
#
# Parts:
#   A. Per-draw time (20 years) and peak memory per grid x parameter set, cutoff auto (full domain) and k = 10
#      (grid cropped to the bounding box of the entry cells + R_c).
#   B. Sensitivity of Norway time-to-event to the cutoff factor k (10 km grid).
#   C. FFT vs V21 across grids: D. virgifera inputs resampled to a 25 km EPSG:3035 grid, invaded area at t = 18
#      vs the V21 baseline (2,010,382 km2). Not an equality test (lon/lat is not shift-invariant).
#
# Usage (package root): Rscript data-raw/benchmark_dispk.R
# Output: data-raw/benchmark_dispk_results.csv (A), data-raw/benchmark_dispk_cutoff.csv (B), printed summary (C).

suppressPackageStartupMessages({ library(terra); devtools::load_all(quiet = TRUE) })

YEARS     <- 20
GRIDS_KM  <- c(20, 10, 5)
EXT_LL    <- terra::ext(-11, 45, 45, 72)          # Europe north of 45 N (lon / lat)
NOR_EXT   <- terra::ext(4, 31.5, 57.5, 71.5)
PARAMS <- list(dvv_like = list(u = 80, p = 5, lmax = 40),
               local    = list(u = 15, p = 50, lmax = 5),
               fat_tail = list(u = 20, p = 2, lmax = 10))
ENTRY_LL <- cbind(c(9.5, 10, 10.5, 9.8, 10.2), c(53.4, 53.6, 53.5, 53.2, 53.8))   # cluster near Hamburg

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
v_eu  <- terra::project(terra::crop(terra::vect(world), EXT_LL), "EPSG:3035")
v_nor <- terra::project(terra::crop(terra::vect(world[which(world$sovereignt == "Norway"), ]), NOR_EXT), "EPSG:3035")

land_grid <- function(res_km) {
  m <- res_km * 1000
  e <- terra::ext(v_eu)
  tmpl <- terra::rast(terra::ext(floor(e$xmin / m) * m, ceiling(e$xmax / m) * m, floor(e$ymin / m) * m,
                                 ceiling(e$ymax / m) * m), resolution = m, crs = "EPSG:3035")
  terra::rasterize(v_eu, tmpl, field = 1, touches = TRUE)
}

entry_xy <- terra::geom(terra::project(terra::vect(ENTRY_LL, crs = "EPSG:4326"), "EPSG:3035"))[, c("x", "y")]

run_one <- function(pra, init, pr, cutoff, k = 10) {
  gc(reset = TRUE)
  tm <- system.time(x <- suppressWarnings(spread_dispk(pra, init, lmax = pr$lmax, u = pr$u, p = pr$p, t = c(5, 10, YEARS),
                                                        cutoff = cutoff, k = k)))
  mem <- sum(gc()[, 6])                                             # max used, Mb (Ncells + Vcells)
  list(x = x, sec = tm[["elapsed"]], mem = mem)
}

# ---- A. per-draw time ----
rows <- list()
for (g in GRIDS_KM) {
  land <- land_grid(g)
  pra <- prep_suitability(land * 0.5, 0.1)
  init <- init_density(pra, entry_xy, 10, to = "suitable")
  nor <- terra::rasterize(v_nor, land, field = 1, touches = TRUE)
  for (nm in names(PARAMS)) {
    pr <- PARAMS[[nm]]
    for (cf in c("auto", "k")) {
      pra_run <- pra; init_run <- init; nor_run <- nor
      if (cf == "k") {                                              # crop to entry bbox + R_c = 10 u
        e <- terra::ext(min(entry_xy[, 1]), max(entry_xy[, 1]), min(entry_xy[, 2]), max(entry_xy[, 2])) + 10 * pr$u * 1000
        e <- terra::intersect(e, terra::ext(land))
        pra_run <- pra; pra_run$rast <- terra::crop(pra$rast, e)
        init_run <- terra::crop(init, e); nor_run <- terra::crop(nor, e)
      }
      res <- run_one(pra_run, init_run, pr, cf)
      x <- res$x
      ok <- identical(x$status, "ok")
      te <- if (ok) time_to_event(x$years_to_reach, pra_run, region = nor_run, pct = c(10, 50), horizon = YEARS) else NULL
      rows[[length(rows) + 1]] <- data.frame(
        grid_km = g, params = nm, u = pr$u, p = pr$p, lmax = pr$lmax, cutoff = cf,
        cells = terra::ncell(pra_run$rast), fft = if (ok) paste(x$params$fft_dim, collapse = "x") else NA,
        R_c_km = if (ok) round(x$params$R_c_km) else NA, kernel_weight_sum = if (ok) round(x$params$kernel_weight_sum, 4) else NA,
        sec_per_draw = round(res$sec, 2), sec_per_year = round(res$sec / YEARS, 3), peak_mem_mb = round(res$mem),
        status = if (ok) "ok" else x$status,
        abundant_km2_t20 = if (ok) round(x$summary$area_abundant_km2[3]) else NA,
        norway_first_reach = if (ok) te$first_year else NA, norway_reach_10pct = if (ok) te$year_pct_10 else NA)
      cat(sprintf("%2d km %-8s %-4s cells %7d  %6.1f s/draw  %5.0f Mb  %s\n", g, nm, cf, terra::ncell(pra_run$rast),
                  res$sec, res$mem, rows[[length(rows)]]$status))
    }
  }
}
A <- do.call(rbind, rows)
write.csv(A, file.path("data-raw", "benchmark_dispk_results.csv"), row.names = FALSE)

# ---- B. cutoff sensitivity (10 km, two parameter sets) ----
land <- land_grid(10)
pra <- prep_suitability(land * 0.5, 0.1)
init <- init_density(pra, entry_xy, 10, to = "suitable")
nor <- terra::rasterize(v_nor, land, field = 1, touches = TRUE)
B <- do.call(rbind, lapply(c("dvv_like", "fat_tail"), function(nm) {
  pr <- PARAMS[[nm]]
  do.call(rbind, lapply(list(3, 5, 10, 20, "auto"), function(kk) {
    x <- suppressWarnings(spread_dispk(pra, init, lmax = pr$lmax, u = pr$u, p = pr$p, t = YEARS,
                                       cutoff = if (identical(kk, "auto")) "auto" else "k", k = if (is.numeric(kk)) kk else 10))
    te <- time_to_event(x$years_to_reach, pra, region = nor, pct = c(10, 50), horizon = YEARS)
    ta <- time_to_event(x$years_to_abundant, pra, region = nor, pct = 10, horizon = YEARS)
    data.frame(params = nm, k = as.character(kk), R_c_km = round(x$params$R_c_km), tail_mass = signif(x$params$tail_mass_beyond_cutoff, 3),
               lost_kernel_share_yr1 = signif(x$mass$lost_kernel[1] / x$mass$dispersing[1], 3),
               nor_first_reach = te$first_year, nor_reach_10pct = te$year_pct_10, nor_reach_50pct = te$year_pct_50,
               nor_first_abundant = ta$first_year, sec = round(x$params$seconds, 1))
  }))
}))
write.csv(B, file.path("data-raw", "benchmark_dispk_cutoff.csv"), row.names = FALSE)

# ---- C. FFT vs V21 across grids (D. virgifera, t = 18) ----
b <- readRDS(file.path("tests", "testthat", "fixtures", "v21_dvv_baseline.rds"))
p_ll <- prep_climex(b$tab, rescale_gi = FALSE, area_method = "v21")
# project the raster itself (projecting only the corners of the lon/lat extent cuts off the LAEA footprint)
ei <- terra::project(p_ll$rast[["ei"]], "EPSG:3035", res = 25000, method = "near")
tmpl <- terra::rast(ei)
gi <- terra::project(p_ll$rast[["gi"]], tmpl, method = "near")
p_ea <- new_pra(ei, gi, area_method = "terra", settings = list(source = "climex resampled", gi_rescale = "domain"))
pres <- b$info$params$presence
init_ea <- init_density(p_ea, as.matrix(pres[, 1:2]), pres[, 3], crs = "EPSG:4326")
tm <- system.time(x <- spread_dispk(p_ea, init_ea, lmax = 40, u = 80, p = 5, t = c(16, 17, 18), cutoff = "auto"))
v21_area <- sum(b$grid$area_km2[b$dispk[, 18] > 1])
cat("\n== C. FFT (25 km EPSG:3035) vs V21 (0.5 deg lon/lat), D. virgifera t = 16..18 ==\n")
cat(sprintf("suitable area: V21 %.0f km2, resampled %.0f km2\n", sum(b$grid$area_km2[b$grid$EI > 0]),
            x$summary$area_suitable_km2[1]))
print(data.frame(t = 16:18, v21_abundant_km2 = sapply(16:18, function(tt) round(sum(b$grid$area_km2[b$dispk[, tt] > 1]))),
                 fft_abundant_km2 = round(x$summary$area_abundant_km2)), row.names = FALSE)
cat(sprintf("t = 18 difference: %+.1f %% of the V21 area; FFT %.1f s (V21 590 s, reference port 196 s)\n",
            100 * (x$summary$area_abundant_km2[3] / v21_area - 1), tm[["elapsed"]]))
# same grid: FFT vs the reference implementation
ref <- spread_dispk_reference(p_ea, init_ea, lmax = 40, u = 80, p = 5, t = c(16, 17, 18))
d1 <- terra::values(ref$density); d2 <- terra::values(x$density); s <- !is.na(d1) & d1 > 0
cat(sprintf("same 25 km grid, FFT vs reference: max rel diff %.1e, abundant km2 t18 %.0f vs %.0f\n",
            max(abs(d1[s] - d2[s]) / d1[s]), x$summary$area_abundant_km2[3], ref$summary$area_abundant_km2[3]))
# N0 is percent of the CELL's K: keeping the founder number (100) on the smaller entry cell
a_entry <- terra::values(p_ea$rast[["area_km2"]])[which(terra::values(init_ea) > 0)]
init_f <- init_density(p_ea, which(terra::values(init_ea) > 0), n0_from(100, k_carrying(a_entry, 0.2, 200e6)))
y <- spread_dispk(p_ea, init_f, lmax = 40, u = 80, p = 5, t = c(16, 17, 18))
cat(sprintf("founders kept (entry cell %.0f vs %.0f km2): abundant km2 %s; t18 %+.1f %% of V21\n",
            a_entry, b$grid$area_km2[which(b$presence > 0)], paste(round(y$summary$area_abundant_km2), collapse = " / "),
            100 * (y$summary$area_abundant_km2[3] / v21_area - 1)))

cat("\n== A. per-draw time ==\n"); print(A, row.names = FALSE)
cat("\n== B. cutoff sensitivity (10 km) ==\n"); print(B, row.names = FALSE)
