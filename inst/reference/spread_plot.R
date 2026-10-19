#
# Author ================================================================================
#
#  Daniel Flo, PhD
#  Panel Coordinator, VKM Panel on Plant Health
#  Norwegian Scientific Committee for Food and Environment (VKM)
#  Email:   daniel.flo@vkm.no
#
# Last updated: 2026-10-19 13:29
# Version: 43
# Sibling: bsspread_europe.R (model-based: can it reach Norway naturally, and how fast).
#
# Intro ================================================================================
#
# Purpose: from European GBIF records, how fast is a species actually coming towards Norway, and when would it arrive
# if the trend continues? Empirical counterpart to the distance screen and simulations in bsspread_europe.R.
#
# Per species (same species list and skip rule as bsspread_europe.R):
#   1. GBIF records in Europe since YEAR_FROM (exact species match), cleaned; only records on land in Europe north of
#      45 N (the bsspread_europe.R grid) count.
#   2. Established records, resolved in time: a record counts as a source from the year its 20 km cell plus neighbours
#      has records in >= 2 different years (the bsspread_europe.R rule). Species without established records: skipped
#      (one table row, no figure).
#   3. Approach rate: distance from Norway to the nearest established record, year by year; approach speed = minus its
#      slope; years to Norway = current distance / approach speed (straight line; bsspread_europe.R uses the route
#      over suitable land, so its screen years are listed alongside for comparison).
#   4. Front speed: distance of each 20 km cell from the origin (earliest established records) regressed on the year
#      the cell became established (Gilbert & Liebhold 2010: most reliable under sparse / clustered sampling); boundary
#      displacement (yearly largest distance, cf. Tobin et al. 2007) as a secondary column.
#   5. Sampling effort: GBIF records per year of the species' ORDER in Europe. 'Effort-corrected' = each replicate keeps
#      records of year y with probability effort_ref / effort_y (all years thinned to the lowest effort in the window);
#      'all records' = no thinning. Both from REPLICATES bootstrap replicates (median, 95% interval).
#   The window starts when the order's effort first reaches EFFORT_MIN_SHARE of its maximum (earlier years have too
#   little recording to show a front) and when the species has established records.
# Outputs: per species in its species folder, as bsspread_europe.R: <Species>/<EPPO>/spread_plot/figures/approach_<EPPO>.png
# and tables/spread_towards_norway_<EPPO>.xlsx (yearly series, established records); species without a species folder
# under 5_outputs/spread_plot/species/<EPPO>/. One table for all species: 5_outputs/spread_plot/spread_towards_norway.xlsx.
# Skipped species get no folder (one row in the table only).
# Approach and front speed include human-assisted spread: they are NOT natural dispersal rates for the trait table.
#
# ========================================================================================

# ---------- CHAPTER 0: CONFIGURATION ----------
SPECIES          <- "ALL"        # "ALL" = every species with an SDM raster (as bsspread_europe.R), or e.g. c("HESOCA", "PAYSAR")
YEAR_FROM        <- 1970         # GBIF records from this year (the approach needs a long series)
LAST_YEAR        <- as.integer(format(Sys.Date(), "%Y")) - 1   # the current year is incomplete
DOWNLOAD_LIMIT   <- 100000       # GBIF records per species
MAX_UNC_M        <- 10000        # drop records with coordinate uncertainty above this (m)
EST_CELL_KM      <- 20           # 'established' rule, as bsspread_europe.R
EST_MIN_YEARS    <- 2
GRID_KM          <- 5            # land grid north of 45 N and distance-to-Norway raster
EFFORT_MIN_SHARE <- 0.1          # window starts when the order's yearly GBIF records reach this share of their maximum
MIN_WINDOW_YEARS <- 10           # fewer years in the window: estimates flagged as unreliable
MIN_EST_THINNED  <- 5            # effort-corrected replicates with fewer established records: use 'all records' instead
REPLICATES       <- 200
SEED             <- 1

DATA_DIR      <- "C:/Users/dafl/OneDrive - Folkehelseinstituttet/Prosjektdata - Dokumenter/VKM Data/26.08.2024_lopende_oppdrag_plantehelse/data"
SPECIES_BASES <- c("C:/Users/dafl/OneDrive - Folkehelseinstituttet/VKM Data/26.08.2024_lopende_oppdrag_plantehelse/Species",
                   "C:/Users/dafl/OneDrive - Folkehelseinstituttet/Prosjektdata - Dokumenter/VKM Data/26.08.2024_lopende_oppdrag_plantehelse/Species")
MASTER_FILE   <- file.path(DATA_DIR, "2_processed_data", "master_species.xlsx")
SCREEN_FILE   <- file.path(DATA_DIR, "5_outputs", "spread_europe", "europe_spread_screen.xlsx")   # optional, for comparison
OUT_DIR       <- file.path(DATA_DIR, "5_outputs", "spread_plot")
SCRIPT_DIR    <- file.path(DATA_DIR, "4_scripts")

suppressPackageStartupMessages({ library(terra); library(sf); library(ggplot2); library(data.table) })
source(file.path(SCRIPT_DIR, "plot_scales.R"))
source(file.path(SCRIPT_DIR, "spread_helpers.R"))
source(file.path(SCRIPT_DIR, "spread_europe_helpers.R"))
source(file.path(SCRIPT_DIR, "approach_helpers.R"))

for (d_ in c("gbif_cache", "effort_cache")) dir.create(file.path(OUT_DIR, d_), showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(SCRIPT_DIR, "logs", paste0("spread_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
log_line <- function(...) {
  msg <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), paste0(...))
  cat(msg, "\n", sep = ""); cat(msg, "\n", sep = "", file = log_file, append = TRUE)
}

# ---------- CHAPTER 1: SHARED INPUTS ----------
log_line("spread_plot.R | log ", log_file)
world  <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
land   <- se_land_grid(world, GRID_KM)
nor    <- se_norway_mask(world, land)
dist_r <- ap_distance_raster(nor)
dist_v <- values(dist_r, mat = FALSE); land_v <- values(land, mat = FALSE)
world_3035  <- st_transform(world, AP_CRS)
norway_3035 <- st_transform(world[which(world$sovereignt == "Norway"), ], AP_CRS)
master <- as.data.frame(readxl::read_excel(MASTER_FILE))
screen <- if (file.exists(SCREEN_FILE)) as.data.frame(readxl::read_excel(SCREEN_FILE)) else NULL
species <- if (identical(SPECIES, "ALL")) se_species_with_sdm(SPECIES_BASES) else SPECIES
log_line("grid: ", GRID_KM, " km land north of 45 N | species: ", length(species), " | years ", YEAR_FROM, "-", LAST_YEAR)

results <- list()
row_df <- function(x) as.data.frame(lapply(x, function(v) if (is.null(v) || length(v) == 0) NA else v),
                                    stringsAsFactors = FALSE, check.names = FALSE)
save_results <- function() {
  cols <- unique(unlist(lapply(results, names)))
  tab <- do.call(rbind, lapply(results, function(x) { x <- x[cols]; names(x) <- cols; row_df(x) }))
  writexl::write_xlsx(list(spread_towards_norway = tab,
                           settings = data.frame(setting = c("YEAR_FROM", "LAST_YEAR", "MAX_UNC_M", "EST_CELL_KM", "EST_MIN_YEARS",
                                                             "GRID_KM", "EFFORT_MIN_SHARE", "MIN_WINDOW_YEARS", "MIN_EST_THINNED",
                                                             "REPLICATES", "SEED"),
                                                 value = c(YEAR_FROM, LAST_YEAR, MAX_UNC_M, EST_CELL_KM, EST_MIN_YEARS, GRID_KM,
                                                           EFFORT_MIN_SHARE, MIN_WINDOW_YEARS, MIN_EST_THINNED, REPLICATES, SEED))),
                      file.path(OUT_DIR, "spread_towards_norway.xlsx"))
}
# median (2.5% - 97.5%); infinite values (not approaching) shown as words
fmt <- function(q, d = 1) {
  if (!is.finite(q[1])) return("not approaching")
  hi <- if (is.finite(q[3])) formatC(q[3], format = "f", digits = d) else "never"
  paste0(formatC(q[1], format = "f", digits = d), " (", formatC(q[2], format = "f", digits = d), "-", hi, ")")
}

# ---------- CHAPTER 2: PER SPECIES ----------
for (i in seq_along(species)) {
  eppo <- species[i]; t_sp <- Sys.time()
  row <- tryCatch({
    sp <- master[!is.na(master$eppocode) & master$eppocode == eppo, , drop = FALSE]
    if (nrow(sp) == 0) stop("EPPO code not in the master species list")
    name <- sp$fullname[1]
    skip <- function(why) { log_line("[", i, "/", length(species), "] ", eppo, " ", name, " | ", why)
                            list(eppo = eppo, species = name, result = why) }
    local({
      # 2.1 Records on land north of 45 N, distance to Norway, established year
      g <- se_gbif_key(name)
      if (is.na(g$key)) return(skip(paste0("skipped: no usable GBIF match (", g$note, ")")))
      raw <- se_download(g$key, file.path(OUT_DIR, "gbif_cache", paste0("gbif_europe_", YEAR_FROM, "_", eppo, ".rds")),
                         YEAR_FROM, DOWNLOAD_LIMIT)
      cl <- se_clean(raw, MAX_UNC_M)
      r <- ap_records(cl$data)[year <= LAST_YEAR]
      cell <- cellFromXY(land, as.matrix(r[, .(x, y)]))
      r <- r[!is.na(cell) & !is.na(land_v[cell])]
      if (!nrow(r)) return(skip("skipped: no GBIF records in Europe north of 45 N"))
      r[, km := dist_v[cellFromXY(land, as.matrix(r[, .(x, y)]))]]
      r[, est_year := ap_est_year(r, EST_CELL_KM, EST_MIN_YEARS)]
      if (all(is.na(r$est_year))) return(skip(paste0("skipped: ", nrow(r), " record(s) north of 45 N, none established")))

      # 2.2 Window: effort of the order, first established year
      b <- tryCatch(rgbif::name_backbone(name = name), error = function(e) NULL)
      order_key <- if (!is.null(b) && !is.null(b$orderKey)) b$orderKey else NA
      effort <- if (is.na(order_key)) NULL else
        ap_effort(order_key, YEAR_FROM:LAST_YEAR, file.path(OUT_DIR, "effort_cache", paste0("effort_order_", order_key, ".rds")))
      eff_start <- if (is.null(effort)) YEAR_FROM else min(effort$year[effort$records >= EFFORT_MIN_SHARE * max(effort$records)])
      w0 <- max(eff_start, min(r$est_year, na.rm = TRUE)); years <- w0:LAST_YEAR
      keep <- if (is.null(effort)) NULL else {
        kp <- ap_keep_prob(effort, years)
        all_y <- min(r$year):LAST_YEAR
        setNames(ifelse(all_y < w0, kp[1], kp[match(all_y, years)]), all_y)   # before the window: the first year's rate
      }

      # 2.3 Estimates: full data (point values, series), replicates for the intervals
      full <- ap_estimate(r, years, EST_CELL_KM, EST_MIN_YEARS)
      origin <- ap_origin(r$est_year, r$x, r$y)
      d_now <- full$D[length(full$D)]
      rep_all  <- ap_replicates(r, years, NULL, REPLICATES, SEED, EST_CELL_KM, EST_MIN_YEARS, origin)
      rep_thin <- if (is.null(keep)) NULL else ap_replicates(r, years, keep, REPLICATES, SEED, EST_CELL_KM, EST_MIN_YEARS, origin)
      s_all <- ap_summary(rep_all, d_now)
      s_thin <- if (is.null(rep_thin)) NULL else ap_summary(rep_thin, d_now)
      thin_ok <- !is.null(rep_thin) && stats::median(rep_thin$n_est) >= MIN_EST_THINNED
      main <- if (thin_ok) s_thin else s_all
      basis <- if (thin_ok) "effort-corrected" else if (is.null(rep_thin)) "all records (no effort data)" else
        "all records (too few records after effort thinning)"
      flags <- c(if (length(years) < MIN_WINDOW_YEARS) paste0("short window (", length(years), " yr)"),
                 if (main$share_approaching < 0.9 && main$share_approaching > 0.1) "approach uncertain (direction varies across replicates)")
      status <- if (d_now == 0) "already established in Norway" else if (main$share_approaching <= 0.1) "not approaching Norway" else
        "approaching Norway"
      scr <- if (!is.null(screen) && eppo %in% screen$eppo) screen[screen$eppo == eppo, , drop = FALSE][1, ] else NULL
      log_line("[", i, "/", length(species), "] ", eppo, " ", name, " | ", nrow(r), " records, ", sum(!is.na(r$est_year)),
               " established | window ", w0, "-", LAST_YEAR, " | now ", round(d_now), " km from Norway | approach ",
               fmt(main$approach), " km/yr | years to Norway ", fmt(main$years, 0), " | front ", fmt(main$front), " km/yr | ", basis)

      # 2.4 Per-species outputs in the species folder, as bsspread_europe.R: <Species>/<EPPO>/spread_plot/ (figures,
      # tables). Species without a species folder (no SDM): 5_outputs/spread_plot/species/<EPPO>/.
      sdm <- se_find_sdm(eppo, SPECIES_BASES)
      sp_dir <- if (!is.null(sdm)) file.path(sdm$base, eppo, "spread_plot") else file.path(OUT_DIR, "species", eppo)
      for (d_ in c("figures", "tables")) dir.create(file.path(sp_dir, d_), showWarnings = FALSE, recursive = TRUE)
      band <- function(m, f) if (is.null(m) || !nrow(m)) NA else apply(m, 2, f)
      rep_main <- if (thin_ok) rep_thin else rep_all
      writexl::write_xlsx(list(
        series = data.frame(year = years, km_to_norway_all_records = round(full$D, 1),
                            km_to_norway_median = round(band(rep_main$D, function(z) stats::median(z, na.rm = TRUE)), 1),
                            km_to_norway_lo = round(band(rep_main$D, function(z) stats::quantile(z, 0.025, na.rm = TRUE)), 1),
                            km_to_norway_hi = round(band(rep_main$D, function(z) stats::quantile(z, 0.975, na.rm = TRUE)), 1),
                            front_boundary_km_from_origin = round(full$F, 1),
                            effort_keep_probability = if (is.null(keep)) NA else round(unname(keep[as.character(years)]), 4)),
        established_records = r[!is.na(est_year), .(x_3035 = round(x), y_3035 = round(y), year, est_year, km_to_norway = round(km, 1))]),
        file.path(sp_dir, "tables", paste0("spread_towards_norway_", eppo, ".xlsx")))

      # 2.5 Figure
      cap <- paste0("GBIF ", YEAR_FROM, "-", LAST_YEAR, ", cleaned, land north of 45 N; established = records in >= ", EST_MIN_YEARS,
                    " years within ~", EST_CELL_KM, " km. Straight-line distances. Speeds include human-assisted spread.\n",
                    "Effort correction: records thinned to the lowest yearly GBIF effort of the order in ", w0, "-", LAST_YEAR,
                    "; ", REPLICATES, " bootstrap replicates; median (95% interval).")
      sub <- paste0(status, " | now ", round(d_now), " km away | approach ", fmt(main$approach), " km/yr | years to Norway ",
                    fmt(main$years, 0), "\nfront speed ", fmt(main$front), " km/yr | ", basis,
                    if (length(flags)) paste0(" | ", paste(flags, collapse = "; ")))
      ggsave(file.path(sp_dir, "figures", paste0("approach_", eppo, ".png")),
             ap_plot(r[!is.na(est_year)], years, full$D, if (thin_ok) rep_thin else rep_all, main$approach[1],
                     paste0(name, ": spread towards Norway"), sub, cap, world_3035, norway_3035),
             width = 14, height = 7.5, dpi = 200, bg = "white")

      list(eppo = eppo, species = name, result = status, basis = basis, flags = paste(flags, collapse = "; "),
           gbif_match = g$note, records = nrow(r), established_records = sum(!is.na(r$est_year)), order_key = order_key,
           window_start = w0, window_end = LAST_YEAR, km_to_norway_now = round(d_now, 1),
           approach_km_yr = round(main$approach[1], 2), approach_lo = round(main$approach[2], 2), approach_hi = round(main$approach[3], 2),
           share_replicates_approaching = round(main$share_approaching, 2),
           years_to_norway = round(main$years[1], 1), years_lo = round(main$years[2], 1), years_hi = round(main$years[3], 1),
           front_km_yr = round(main$front[1], 2), front_lo = round(main$front[2], 2), front_hi = round(main$front[3], 2),
           front_boundary_km_yr = round(main$front_boundary[1], 2), front_boundary_lo = round(main$front_boundary[2], 2),
           front_boundary_hi = round(main$front_boundary[3], 2),
           approach_km_yr_all_records = round(s_all$approach[1], 2), front_km_yr_all_records = round(s_all$front[1], 2),
           screen_route_km = if (!is.null(scr)) scr$nearest_route_km else NA,
           screen_years_natural_most_likely = if (!is.null(scr)) scr$years_natural_screen else NA,
           screen_class = if (!is.null(scr)) scr$class else NA,
           output_folder = sp_dir,
           minutes = round(as.numeric(difftime(Sys.time(), t_sp, units = "mins")), 2), run_date = format(Sys.time(), "%Y-%m-%d %H:%M"))
    })
  }, error = function(e) {
    log_line("[", i, "/", length(species), "] ", eppo, " FAILED: ", conditionMessage(e))
    list(eppo = eppo, result = paste("FAILED:", conditionMessage(e)))
  })
  results[[eppo]] <- row
  save_results()
}

# ---------- CHAPTER 3: SUMMARY ----------
res <- vapply(results, function(x) sub(":.*", "", x$result), character(1))
log_line("FINISHED | ", length(results), " species | ", file.path(OUT_DIR, "spread_towards_norway.xlsx"))
for (k in names(sort(table(res), decreasing = TRUE))) log_line("  ", sum(res == k), " x ", k)
