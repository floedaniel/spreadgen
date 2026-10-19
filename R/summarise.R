# Tidy summaries of Monte Carlo runs: quantiles and mean per metric; formatted intervals.

#' Summarise Monte Carlo runs
#'
#' @param x An `sg_mc` object or a list of them.
#' @param probs Quantile probabilities (default 0.05, 0.25, 0.5, 0.75, 0.95).
#' @return List with `metrics` (one row per identifier x model x t x level x metric) and `events` (one row per
#'   identifier x model x level x event), each with `n`, `mean`, `q05` ... `q95`; events also `share_never`
#'   (share of successful draws without the event within the horizon; quantiles may be `Inf`).
#' @export
summarise_spread <- function(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  if (inherits(x, "sg_mc")) x <- list(x)
  stopifnot(all(vapply(x, inherits, logical(1), "sg_mc")))
  qn <- paste0("q", formatC(round(100 * probs), width = 2, flag = "0"))
  summ <- function(v) {
    v <- v[!is.na(v)]
    q <- if (length(v)) stats::quantile(v, probs, names = FALSE, type = 7) else rep(NA_real_, length(probs))
    q[is.nan(q)] <- Inf
    c(n = length(v), mean = if (length(v)) mean(v) else NA_real_, stats::setNames(q, qn))
  }
  long <- function(df, keys, vals, extra = NULL) {
    if (is.null(df) || !nrow(df)) return(NULL)
    grp <- interaction(df[keys], drop = TRUE, lex.order = TRUE)
    do.call(rbind, lapply(split(seq_len(nrow(df)), grp), function(ix) {
      head <- df[ix[1], keys, drop = FALSE]
      do.call(rbind, lapply(vals, function(v) {
        z <- df[[v]][ix]
        out <- cbind(head, metric = v, as.data.frame(as.list(summ(z))), row.names = NULL)
        if (!is.null(extra)) out <- cbind(out, extra(z))
        out
      }))
    }))
  }
  id_cols <- function(df) setdiff(names(df), c("draw", "t", "level", "n_cells", "area_km2", "pct_suitable",
                                               "target_area_km2", "target_pct_suitable", "max_distance_km",
                                               grep("year", names(df), value = TRUE)))
  metrics <- do.call(rbind, lapply(x, function(m) {
    df <- m$metrics
    if (is.null(df)) return(NULL)
    df$level[is.na(df$level)] <- -1                 # RR: no density level (arrival); NA cannot be a split key
    out <- long(df, c(id_cols(df), "t", "level"),
                c("n_cells", "area_km2", "pct_suitable", "target_area_km2", "target_pct_suitable", "max_distance_km"))
    out$level[out$level == -1] <- NA
    out
  }))
  events <- do.call(rbind, lapply(x, function(m) {
    df <- m$events
    if (is.null(df)) return(NULL)
    df$level[is.na(df$level)] <- -1
    out <- long(df, c(id_cols(df), "level"), grep("year", names(df), value = TRUE),
                extra = function(z) data.frame(share_never = mean(is.infinite(z))))
    out$level[out$level == -1] <- NA
    out
  }))
  rownames(metrics) <- NULL
  rownames(events) <- NULL
  list(metrics = metrics, events = events)
}

#' Format a sample as "median (2.5%-97.5%)"
#'
#' @param x Numeric sample (may contain `Inf` = no event within the horizon).
#' @param digits Decimal places.
#' @param never Word for infinite values.
#' @param probs Lower and upper probability.
#' @return Character string, e.g. `"12 (8-never within horizon)"`; `NA` when `x` is empty.
#' @export
fmt_interval <- function(x, digits = 0, never = "never within horizon", probs = c(0.025, 0.975)) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  q <- stats::quantile(x, c(0.5, probs), names = FALSE, type = 7)
  q[is.nan(q)] <- Inf
  f <- function(v) if (is.finite(v)) formatC(v, format = "f", digits = digits) else never
  if (!is.finite(q[1])) return(never)
  paste0(f(q[1]), " (", f(q[2]), "-", f(q[3]), ")")
}
