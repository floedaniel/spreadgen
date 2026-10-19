# M0 baseline: run the original Prog_SpreadModule_V21-web.R (Robinet et al. 2012; DANS doi:10.17026/DANS-2XV-Z3DX, CC0)
# under modern R on the D. virgifera virgifera case-study inputs, and save the results as regression fixtures.
#
# The original file is sourced unchanged. Minimal patches are shims defined around it:
#   1. xyValues() (removed from 'raster') -> raster::extract() (same argument order: raster, lon/lat matrix).
#   2. plotkernel() (called every year inside dispk()) is replaced, after sourcing, by a recorder that stores the yearly
#      density vector; no plotting inside the computation loop.
#   3. All other base-graphics output goes to a null device.
#   4. 'raster' is attached before sourcing (the original expected the user to load sp/raster/rgdal by hand; rgdal is
#      not needed for the functions run here).
# Interactive globals (climexcsv, elevmax, habitatfile, habitatformat) are set as in the case study's "useful commands.r".
#
# Parameters (tutorial / paper, D. virgifera): RR = 80 km/yr, entry Belgrade (20.30 E, 44.82 N), lmax = 40,
# N0 = 100 * 100 / 6.3e10 % (SLG), DK: u = 80 km, p = 5, presence.txt (lon, lat, density %), t = 1..18.
#
# Usage (from the package root, after data-raw/download_original.R): Rscript data-raw/run_v21_baseline.R
# Output: tests/testthat/fixtures/v21_dvv_baseline.rds

suppressPackageStartupMessages(library(raster))

T_MAX  <- 18
RR     <- 80
COORD  <- c(20.30, 44.82)
LMAX   <- 40
N0_SLG <- 100 * 100 / 6.3e10
U      <- 80
P_DF   <- 5

src_dir <- normalizePath(file.path("inst", "original", "Code for case studies", "1-Diabrotica virgifera virgifera"))
out     <- normalizePath(file.path("tests", "testthat", "fixtures"), mustWork = FALSE)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# Work in a temporary copy: the original reads and writes relative to the working directory
work <- file.path(tempdir(), "v21_dvv")
dir.create(work, showWarnings = FALSE)
invisible(file.copy(list.files(src_dir, full.names = TRUE), work, overwrite = TRUE))
old_wd <- setwd(work)
on.exit(setwd(old_wd), add = TRUE)

# Shim 1 and globals expected by the original
xyValues      <- function(object, xy, ...) raster::extract(object, xy)
climexcsv     <- FALSE
elevmax       <- FALSE
habitatfile   <- TRUE
habitatformat <- ".tif"

grDevices::pdf(NULL)                                   # shim 3
t_load <- system.time(source("Prog_SpreadModule_V21-web.R", echo = FALSE, print.eval = FALSE))

# Shim 2: record dispk() yearly states instead of plotting them
dk_years <- list()
plotkernel <- function(vectpop0, t = NULL, ...) dk_years[[t]] <<- vectpop0

# ---- Model C: SLG, t = 1..T_MAX (closed form, so one call per t) ----
slg_runs <- lapply(seq_len(T_MAX), function(t) slg(N0 = N0_SLG, lmax = LMAX, t = t))
slg_mat  <- sapply(slg_runs, `[[`, "slg")
slg_sum  <- data.frame(model = "SLG", t = seq_len(T_MAX),
                       nRA = sapply(slg_runs, `[[`, "nRA"), ntot = sapply(slg_runs, `[[`, "ntot"),
                       sum0 = sapply(slg_runs, `[[`, "sum0"), sum25 = sapply(slg_runs, `[[`, "sum25"),
                       sum50 = sapply(slg_runs, `[[`, "sum50"), sum75 = sapply(slg_runs, `[[`, "sum75"),
                       sum100 = sapply(slg_runs, `[[`, "sum100"))

# ---- Model B: radial, t = 1..T_MAX; invaded mask on the CLIMEX point list (tab) ----
key <- function(lat, lon) paste(round(lat, 3), round(lon, 3))
tab_key <- key(tab[, 1], tab[, 2])
t_rad <- system.time(
  rad_runs <- lapply(seq_len(T_MAX), function(t) radial(RR = RR, t = t, coord = COORD))
)
rad_mat <- sapply(rad_runs, function(r) tab_key %in% key(r$radial[, 1], r$radial[, 2]))
rad_sum <- data.frame(model = "RR", t = seq_len(T_MAX),
                      nRA = sapply(rad_runs, `[[`, "nRA"), ninv = sapply(rad_runs, `[[`, "ninv"),
                      pniche = sapply(rad_runs, `[[`, "pniche"))

# ---- Model D: deterministic dispersal kernel, one run to T_MAX, yearly states recorded ----
t_dk <- system.time(dk <- dispk(N0 = NULL, t = T_MAX, lmax = LMAX, p = P_DF, u = U, presencefile = TRUE))
dk_mat <- do.call(cbind, dk_years)
dk_sum <- data.frame(model = "DK", t = seq_len(T_MAX), nRA = sum(!vectexpr0),
                     ninv = colSums(dk_mat > 1), pniche = 100 * colSums(dk_mat > 1) / sum(!vectexpr0),
                     total_mass = colSums(dk_mat * area / 100))   # sum(density% / 100 * cell area km2): occupied-K km2
stopifnot(isTRUE(all.equal(dk_mat[, T_MAX], dk$dispk)), dk$ninv == dk_sum$ninv[T_MAX])

grDevices::dev.off()

baseline <- list(
  info = list(source = "Prog_SpreadModule_V21-web.R (DANS id 23542), case study 1-Diabrotica virgifera virgifera",
              doi = "10.17026/DANS-2XV-Z3DX", R = R.version.string, raster = as.character(packageVersion("raster")),
              run_date = format(Sys.time(), "%Y-%m-%d %H:%M"),
              params = list(T_MAX = T_MAX, RR = RR, coord = COORD, lmax = LMAX, N0_slg = N0_SLG, u = U, p = P_DF,
                            presence = read.table("presence.txt")),
              seconds = c(load = t_load[["elapsed"]], radial = t_rad[["elapsed"]], dispk = t_dk[["elapsed"]])),
  # CLIMEX point list after habitat masking (lat, lon, EI, GI rescaled), used by SLG / radial
  tab = data.frame(lat = tab[, 1], lon = tab[, 2], EI = EI, GI = GI, na_any = exprNATOT),
  # full rectangular grid used by dispk (lon, lat), with EI / GI / area / sea flag
  grid = data.frame(lon = vectxy[, 1], lat = vectxy[, 2], EI = vectEI, GI = vectGI, area_km2 = area, sea = vectsea),
  slg = slg_mat, radial = rad_mat, dispk = dk_mat, presence = dk$presence,
  summary = list(slg = slg_sum, radial = rad_sum, dispk = dk_sum)
)
saveRDS(baseline, file.path(out, "v21_dvv_baseline.rds"), compress = "xz")

cat(sprintf("Grid: %d CLIMEX points, %d cells in dispk grid, nRA = %d, mean cell area %.0f km2\n",
            nrow(tab), nrow(vectxy), sum(EI > 0), mean(area)))
cat(sprintf("Seconds: load %.1f, radial (18 runs) %.1f, dispk (t = %d) %.1f\n",
            t_load[["elapsed"]], t_rad[["elapsed"]], T_MAX, t_dk[["elapsed"]]))
print(rbind(rad_sum[rad_sum$t %in% c(5, 10, 15, 18), c("model", "t", "nRA", "ninv", "pniche")],
            dk_sum[dk_sum$t %in% c(5, 10, 15, 18), c("model", "t", "nRA", "ninv", "pniche")]), row.names = FALSE)
print(slg_sum[slg_sum$t %in% c(5, 10, 15, 18), ], row.names = FALSE)
cat("Saved", file.path(out, "v21_dvv_baseline.rds"), "\n")
