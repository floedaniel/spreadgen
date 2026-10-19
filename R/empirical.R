# Radial-rate priors from an empirical front-speed table (e.g. the spread_plot.R output).

#' Radial-rate (RR) priors from an empirical front-speed table
#'
#' Reads a table with one row per species and the columns `eppo`, `front_km_yr`, `front_lo`, `front_hi` (km/yr;
#' median and 95 % interval) and optionally `flags`, `basis`, `result` (e.g. `spread_towards_norway.xlsx` written by
#' `spread_plot.R`). A row is usable when the three values are finite, the median is positive, `front_lo <= median <=
#' front_hi` and no flag matches `unreliable`. With `nonpositive_lo = "floor"` (default) a lower bound `<= 0` keeps
#' the row usable and sets `floor_needed = TRUE`: the lower bound is then replaced with [rr_floor()]. With
#' `"unusable"` such rows are unusable. Empirical front speeds include human-assisted spread: they are total spread
#' rates, not natural dispersal rates.
#'
#' @param x Path to an .xlsx file (first sheet) or a data frame.
#' @param unreliable Regular expression for flags that make a row unusable (default: short window).
#' @param nonpositive_lo `"floor"` or `"unusable"` (see Details).
#' @return Data frame: `eppo`, `rr_min` (raw lower bound), `rr_mode`, `rr_max`, `usable`, `floor_needed`, `reason`,
#'   `basis`, `flags`, `rr_source = "empirical front speed"`.
#' @export
read_empirical_priors <- function(x, unreliable = "short window", nonpositive_lo = c("floor", "unusable")) {
  nonpositive_lo <- match.arg(nonpositive_lo)
  if (is.character(x)) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("reading .xlsx needs the 'readxl' package")
    x <- as.data.frame(readxl::read_excel(x))
  }
  need <- c("eppo", "front_km_yr", "front_lo", "front_hi")
  miss <- setdiff(need, names(x))
  if (length(miss)) stop("empirical table lacks column(s): ", paste(miss, collapse = ", "))
  num <- function(v) suppressWarnings(as.numeric(v))
  lo <- num(x$front_lo)
  md <- num(x$front_km_yr)
  hi <- num(x$front_hi)
  flags <- if ("flags" %in% names(x)) as.character(x$flags) else rep(NA_character_, nrow(x))
  flags[is.na(flags)] <- ""
  finite <- is.finite(lo) & is.finite(md) & is.finite(hi)
  ordered <- finite & md > 0 & lo <= md & md <= hi
  floor_needed <- ordered & lo <= 0
  if (nonpositive_lo == "unusable") ordered <- ordered & !floor_needed
  flagged <- nzchar(unreliable) & grepl(unreliable, flags)
  reason <- ifelse(!finite, "no finite front speed", ifelse(!ordered, if (nonpositive_lo == "unusable")
                   "front speed interval not positive / ordered" else "front speed median not positive / interval not ordered",
                   ifelse(flagged, paste0("flagged: ", flags), "")))
  data.frame(eppo = as.character(x$eppo), rr_min = lo, rr_mode = md, rr_max = hi, usable = ordered & !flagged,
             floor_needed = floor_needed & !flagged, reason = reason,
             basis = if ("basis" %in% names(x)) as.character(x$basis) else NA_character_,
             flags = flags, rr_source = "empirical front speed", stringsAsFactors = FALSE)
}

#' Floor for a non-positive lower bound of an empirical front speed
#'
#' A lower bound `<= 0` cannot be a PERT minimum on the log scale. It is replaced with the lower bound of the
#' species' natural rate (trait table) when that is positive and below the median, otherwise with `0.1 x median`.
#' Lower bounds `> 0` are kept.
#'
#' @param rr_min,rr_mode Lower bound and median of the empirical front speed (km/yr).
#' @param natural_min Lower bound of the natural rate (trait table; `NA` when there is none).
#' @return Data frame: `rr_min` (floored), `floor_applied` (logical), `floor_source` (`"trait table"`,
#'   `"0.1 x median"` or `NA`).
#' @export
rr_floor <- function(rr_min, rr_mode, natural_min = NA_real_) {
  n <- max(length(rr_min), length(rr_mode), length(natural_min))
  rr_min <- rep_len(as.numeric(rr_min), n)
  rr_mode <- rep_len(as.numeric(rr_mode), n)
  natural_min <- rep_len(as.numeric(natural_min), n)
  apply <- is.finite(rr_min) & rr_min <= 0
  nat_ok <- is.finite(natural_min) & natural_min > 0 & natural_min < rr_mode
  out_min <- ifelse(apply, ifelse(nat_ok, natural_min, 0.1 * rr_mode), rr_min)
  data.frame(rr_min = out_min, floor_applied = apply,
             floor_source = ifelse(apply, ifelse(nat_ok, "trait table", "0.1 x median"), NA_character_),
             stringsAsFactors = FALSE)
}
