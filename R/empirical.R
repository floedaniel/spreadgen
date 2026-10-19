# Radial-rate priors from an empirical front-speed table (e.g. the spread_plot.R output).

#' Radial-rate (RR) priors from an empirical front-speed table
#'
#' Reads a table with one row per species and the columns `eppo`, `front_km_yr`, `front_lo`, `front_hi` (km/yr;
#' median and 95 % interval) and optionally `flags`, `basis`, `result` (e.g. `spread_towards_norway.xlsx` written by
#' `spread_plot.R`). A row is usable when the three values are finite, `front_lo > 0` and no flag matches
#' `unreliable`. Empirical front speeds include human-assisted spread: they are total spread rates, not natural
#' dispersal rates.
#'
#' @param x Path to an .xlsx file (first sheet) or a data frame.
#' @param unreliable Regular expression for flags that make a row unusable (default: short window).
#' @return Data frame: `eppo`, `rr_min`, `rr_mode`, `rr_max`, `usable`, `reason`, `basis`, `flags`,
#'   `rr_source = "empirical front speed"`.
#' @export
read_empirical_priors <- function(x, unreliable = "short window") {
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
  ordered <- finite & lo > 0 & lo <= md & md <= hi
  flagged <- nzchar(unreliable) & grepl(unreliable, flags)
  reason <- ifelse(!finite, "no finite front speed", ifelse(!ordered, "front speed interval not positive / ordered",
                   ifelse(flagged, paste0("flagged: ", flags), "")))
  data.frame(eppo = as.character(x$eppo), rr_min = lo, rr_mode = md, rr_max = hi, usable = ordered & !flagged,
             reason = reason, basis = if ("basis" %in% names(x)) as.character(x$basis) else NA_character_,
             flags = flags, rr_source = "empirical front speed", stringsAsFactors = FALSE)
}
