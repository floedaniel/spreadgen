# ggplot2 figures with viridis scales (scales can be replaced by the caller, e.g. a project's own scale functions).

#' Map of the per-cell probability of a Monte Carlo run
#'
#' @param g An `sg_grid` object.
#' @param prob Numeric vector over all cells (e.g. a column of `sg_mc$prob`).
#' @param title,subtitle,caption Text.
#' @param fill_scale A ggplot2 fill scale (default viridis, limits 0-1).
#' @param borders Optional sf object drawn on top (outlines).
#' @param show Cells to draw: `"suitable"` (default) or `"data"` (all cells with data).
#' @param region `"all"` (default) or `"target"`: draw probabilities in the target cells only, on the target's
#'   bounding box (plus `margin_km`); other land cells in the box are grey.
#' @param margin_km Margin around the target's bounding box (`region = "target"`).
#' @return A ggplot object.
#' @export
plot_probability <- function(g, prob, title = NULL, subtitle = NULL, caption = NULL,
                             fill_scale = ggplot2::scale_fill_viridis_c("Probability", limits = c(0, 1), option = "viridis"),
                             borders = NULL, show = c("suitable", "data"), region = c("all", "target"), margin_km = 25) {
  show <- match.arg(show)
  region <- match.arg(region)
  keep <- if (show == "suitable") which(!is.na(g$ei) & g$ei == 1) else which(!is.na(g$ei))
  lim <- list(x = g$ext[1:2], y = g$ext[3:4])
  if (region == "target") {
    if (!any(g$target)) stop("region = \"target\" needs target cells")
    keep <- keep[g$target[keep]]
    txy <- grid_xy(g, which(g$target))
    m <- margin_km * 1000 * (if (isTRUE(g$lonlat)) 1 / 111320 else 1)
    lim <- list(x = range(txy[, 1]) + c(-1, 1) * (m + g$res[1] / 2), y = range(txy[, 2]) + c(-1, 1) * (m + g$res[2] / 2))
  }
  xy <- grid_xy(g, keep)
  df <- data.frame(x = xy[, 1], y = xy[, 2], p = prob[keep])
  bg <- grid_xy(g, which(!is.na(g$ei)))
  inside <- bg[, 1] >= lim$x[1] & bg[, 1] <= lim$x[2] & bg[, 2] >= lim$y[1] & bg[, 2] <= lim$y[2]
  bg <- bg[inside, , drop = FALSE]
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = data.frame(x = bg[, 1], y = bg[, 2]), ggplot2::aes(.data$x, .data$y), fill = "grey90") +
    ggplot2::geom_raster(data = df, ggplot2::aes(.data$x, .data$y, fill = .data$p)) +
    fill_scale +
    ggplot2::coord_equal(xlim = lim$x, ylim = lim$y, expand = FALSE) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_blank())
  if (!is.null(borders)) {
    p <- p + ggplot2::geom_sf(data = borders, fill = NA, colour = "grey30", linewidth = 0.2, inherit.aes = FALSE) +
      ggplot2::coord_sf(crs = g$crs, xlim = lim$x, ylim = lim$y, expand = FALSE)
  }
  p
}

#' Cumulative distribution of a time-to-event over draws
#'
#' @param events Data frame of events (e.g. `rbind` of `sg_mc$events`), with a grouping column.
#' @param event Event column (default `"target_first_year"`).
#' @param group Column that colours the curves (e.g. `"label"`).
#' @param horizon Last simulated year (draws without the event stay below 1).
#' @param title,subtitle,caption Text.
#' @param colour_scale A ggplot2 colour scale (default viridis, discrete).
#' @return A ggplot object: share of draws in which the event has happened, by year.
#' @export
plot_events <- function(events, event = "target_first_year", group = "label", horizon = NULL, title = NULL,
                        subtitle = NULL, caption = NULL,
                        colour_scale = ggplot2::scale_colour_viridis_d(NULL, option = "viridis", end = 0.9)) {
  if (is.null(horizon)) horizon <- max(events[[event]][is.finite(events[[event]])], 1, na.rm = TRUE)
  yrs <- 0:horizon
  df <- do.call(rbind, lapply(split(events, events[[group]]), function(d) {
    v <- d[[event]]
    data.frame(group = d[[group]][1], year = yrs, share = vapply(yrs, function(y) mean(v <= y, na.rm = TRUE), numeric(1)))
  }))
  ggplot2::ggplot(df, ggplot2::aes(.data$year, .data$share, colour = .data$group)) +
    ggplot2::geom_step(linewidth = 0.9) +
    colour_scale +
    ggplot2::scale_y_continuous("Share of draws", limits = c(0, 1)) +
    ggplot2::labs(x = "Years", title = title, subtitle = subtitle, caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

#' Metric over time: median and quantile bands over draws
#'
#' @param metrics Data frame of metrics (e.g. `rbind` of `sg_mc$metrics`) with a grouping column.
#' @param metric Metric column (default `"target_pct_suitable"`).
#' @param group Column for colour / fill.
#' @param ylab,title,subtitle,caption Text.
#' @param colour_scale,fill_scale ggplot2 scales (default viridis, discrete).
#' @return A ggplot object (median line, 25-75 % and 5-95 % bands).
#' @export
plot_metric <- function(metrics, metric = "target_pct_suitable", group = "label", ylab = metric, title = NULL,
                        subtitle = NULL, caption = NULL,
                        colour_scale = ggplot2::scale_colour_viridis_d(NULL, option = "viridis", end = 0.9),
                        fill_scale = ggplot2::scale_fill_viridis_d(NULL, option = "viridis", end = 0.9)) {
  q <- function(v, p) stats::quantile(v, p, na.rm = TRUE, names = FALSE)
  df <- do.call(rbind, lapply(split(metrics, list(metrics[[group]], metrics$t), drop = TRUE), function(d) {
    v <- d[[metric]]
    data.frame(group = d[[group]][1], t = d$t[1], q05 = q(v, 0.05), q25 = q(v, 0.25), q50 = q(v, 0.5),
               q75 = q(v, 0.75), q95 = q(v, 0.95))
  }))
  ggplot2::ggplot(df, ggplot2::aes(.data$t)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$q05, ymax = .data$q95, fill = .data$group), alpha = 0.15) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$q25, ymax = .data$q75, fill = .data$group), alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = .data$q50, colour = .data$group), linewidth = 0.9) +
    ggplot2::geom_point(ggplot2::aes(y = .data$q50, colour = .data$group), size = 1.6) +
    colour_scale + fill_scale +
    ggplot2::labs(x = "Years", y = ylab, title = title, subtitle = subtitle, caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}
