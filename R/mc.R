# Monte Carlo over PERT-distributed parameters (and sampled entry points), on a plain-data grid.

#' PERT prior for one parameter
#'
#' @param min,mode,max Minimum, most likely and maximum value (`min <= mode <= max`).
#' @param log Sample on the log scale (PERT on log(min), log(mode), log(max)); for rates and multiplication factors.
#' @param source,confidence Free text kept with the prior (e.g. "empirical front speed", "low").
#' @return An object of class `sg_prior`.
#' @export
prior_pert <- function(min, mode, max, log = FALSE, source = NA_character_, confidence = NA_character_) {
  stopifnot(is.numeric(c(min, mode, max)), min <= mode, mode <= max)
  if (log) stopifnot(min > 0)
  structure(list(min = min, mode = mode, max = max, log = log, source = source, confidence = confidence),
            class = "sg_prior")
}

#' PERT prior from a single value (mode) with relative bounds
#'
#' @param x Most likely value.
#' @param rel Relative half-width (default 0.25: min = 0.75 x, max = 1.25 x).
#' @inheritParams prior_pert
#' @export
prior_from_value <- function(x, rel = 0.25, log = TRUE, source = NA_character_, confidence = NA_character_) {
  prior_pert(x * (1 - rel), x, x * (1 + rel), log, source, confidence)
}

#' PERT prior from a range
#'
#' Mode = geometric mean of the bounds when `log = TRUE`, arithmetic mean otherwise.
#' @param lo,hi Bounds.
#' @inheritParams prior_pert
#' @export
prior_from_range <- function(lo, hi, log = TRUE, source = NA_character_, confidence = NA_character_) {
  prior_pert(lo, if (log) sqrt(lo * hi) else (lo + hi) / 2, hi, log, source, confidence)
}

#' @export
print.sg_prior <- function(x, ...) {
  cat(sprintf("<sg_prior> PERT(%g, %g, %g)%s%s\n", x$min, x$mode, x$max, if (x$log) " on log scale" else "",
              if (!is.na(x$source)) paste0(" | ", x$source) else ""))
  invisible(x)
}

#' PERT distribution (beta-PERT, shape lambda = 4)
#'
#' @param n Number of draws.
#' @param p Probabilities.
#' @param min,mode,max Bounds and mode.
#' @param lambda Shape (4 = standard PERT).
#' @return Draws (`rpert`) or quantiles (`qpert`).
#' @export
rpert <- function(n, min, mode, max, lambda = 4) qpert(stats::runif(n), min, mode, max, lambda)

#' @rdname rpert
#' @export
qpert <- function(p, min, mode, max, lambda = 4) {
  if (max == min) return(rep(min, length(p)))
  a <- 1 + lambda * (mode - min) / (max - min)
  b <- 1 + lambda * (max - mode) / (max - min)
  min + stats::qbeta(p, a, b) * (max - min)
}

#' Latin hypercube sample on the unit cube (base R)
#'
#' One stratum per draw and dimension, random position inside the stratum, independent random permutations.
#' @param n Draws.
#' @param k Dimensions.
#' @return Matrix n x k in (0, 1).
#' @export
lhs_unif <- function(n, k) {
  if (k == 0) return(matrix(numeric(0), n, 0))
  vapply(seq_len(k), function(j) (sample.int(n) - stats::runif(n)) / n, numeric(n))
}

#' Draw parameter sets from priors
#'
#' @param priors Named list of `sg_prior` objects.
#' @param n Number of draws.
#' @param method `"lhs"` (Latin hypercube) or `"random"`.
#' @return Data frame, one row per draw.
#' @export
draw_params <- function(priors, n, method = c("lhs", "random")) {
  method <- match.arg(method)
  stopifnot(length(priors) > 0, all(vapply(priors, inherits, logical(1), "sg_prior")))
  U <- if (method == "lhs") lhs_unif(n, length(priors)) else matrix(stats::runif(n * length(priors)), n)
  out <- lapply(seq_along(priors), function(j) {
    pr <- priors[[j]]
    if (pr$log) exp(qpert(U[, j], log(pr$min), log(pr$mode), log(pr$max))) else qpert(U[, j], pr$min, pr$mode, pr$max)
  })
  stats::setNames(as.data.frame(out), names(priors))
}

#' Deterministic seed for one unit of work (e.g. one species)
#'
#' Results then do not depend on the order or the number of parallel workers.
#' @param seed Base seed (integer).
#' @param key Character key (e.g. EPPO code and model); NULL = `seed`.
#' @return Integer seed.
#' @export
seed_from_key <- function(seed, key = NULL) {
  if (is.null(key)) return(as.integer(seed))
  h <- as.numeric(seed) %% 2147483647
  for (cp in utf8ToInt(paste(key, collapse = "|"))) h <- (h * 131 + cp) %% 2147483647
  as.integer(h)
}

#' Monte Carlo run of one spread model
#'
#' Draws parameters from PERT priors (Latin hypercube by default) and, for introduction scenarios, entry cells; runs
#' the model core for every draw and returns tidy per-draw metrics, time-to-event per draw, and the per-cell share of
#' draws at or above the headline level. All randomness is drawn up front from `seed_from_key(seed, key)`; the
#' caller's RNG state is restored.
#'
#' Parameters per model (names of `priors`):
#' * `"SLG"`: `lmax`, `n0` (initial density in every suitable cell, percent of K).
#' * `"RR"`: `RR` (km/yr).
#' * `"DK"`: `RR` (km/yr; the kernel scale is `u = RR`, the paper's rule, unless a `u` prior is given), `p`, `lmax`,
#'   `init` (initial density in the entry cells, percent of K).
#'
#' @param g An `sg_grid` object ([as_grid()]); its `target` vector defines the target area.
#' @param model `"SLG"`, `"RR"` or `"DK"`.
#' @param priors Named list of `sg_prior` objects.
#' @param n_draws Number of draws.
#' @param t Years at which metrics are reported (the horizon is `max(t)`).
#' @param levels Density levels (percent of K) for SLG and DK (headline 1).
#' @param headline Level used for the per-cell probability (must be in `levels`).
#' @param entry Fixed entry: list with `cells` (entry cell numbers) and optionally `dist_km` (distance from every cell
#'   to the nearest entry cell, e.g. computed at the I/O edge; computed with [grid_dist_km()] when absent).
#' @param entry_sampler Introduction scenario: list with `weights` (non-negative, over all cells; 0 = never) and `n`
#'   (entry cells per draw). Used when `entry` is NULL (RR, DK).
#' @param method `"lhs"` or `"random"`.
#' @param seed,key Seed and key for [seed_from_key()].
#' @param pct Shares of the target's suitable area for time-to-event.
#' @param dk_args Further arguments for [core_dispk()] (e.g. `cutoff`, `min_u_cells`).
#' @param id Named list of identifiers copied into every output row (e.g. `list(species = "DIABVI", scenario = "natural")`).
#' @param time_limit Elapsed-time cap in seconds for the whole run (default `Inf`). Checked before every draw, so a
#'   run can exceed it by one draw; draws not started get status `"timeout"` and `settings$timed_out` is `TRUE`.
#' @return An object of class `sg_mc`: `draws` (parameters and `status` per draw), `metrics` (draw x t x level),
#'   `events` (draw x level), `prob` (matrix cells x t: share of successful draws at or above the headline level, or
#'   invaded for RR), `n_ok`, `settings` (including `timed_out`, `elapsed_s`).
#' @export
mc_spread <- function(g, model = c("SLG", "RR", "DK"), priors, n_draws, t = c(5, 10, 20, 50), levels = c(0.1, 1, 10),
                      headline = 1, entry = NULL, entry_sampler = NULL, method = c("lhs", "random"), seed = 1, key = NULL,
                      pct = c(10, 50), dk_args = list(), id = list(), time_limit = Inf) {
  model <- match.arg(model)
  method <- match.arg(method)
  stopifnot(inherits(g, "sg_grid"), n_draws >= 1, headline %in% levels)
  need <- switch(model, SLG = c("lmax", "n0"), RR = "RR", DK = c("RR", "p", "lmax", "init"))
  miss <- setdiff(need, names(priors))
  if (length(miss)) stop("priors for model ", model, " lack: ", paste(miss, collapse = ", "))
  if (model != "SLG" && is.null(entry) && is.null(entry_sampler)) stop("model ", model, " needs entry or entry_sampler")
  priors <- priors[intersect(c(need, "u"), names(priors))]
  s <- seed_from_key(seed, key)
  rnd <- with_seed(s, {
    dr <- draw_params(priors, n_draws, method)
    ent <- if (model != "SLG" && is.null(entry)) {
      w <- entry_sampler$weights
      w[is.na(w) | is.na(g$ei)] <- 0
      if (sum(w) <= 0) stop("entry_sampler weights are all zero")
      lapply(seq_len(n_draws), function(i) sample.int(length(w), entry_sampler$n, replace = FALSE, prob = w))
    } else NULL
    list(draws = dr, entries = ent)
  })
  draws <- rnd$draws
  if (model == "DK" && !"u" %in% names(draws)) draws$u <- draws$RR
  horizon <- max(t)
  suit <- !is.na(g$ei) & g$ei == 1
  g_dom <- g
  g_dom$target <- rep(TRUE, length(g$ei))
  fixed_d <- if (model != "SLG" && !is.null(entry)) {
    if (!is.null(entry$dist_km)) entry$dist_km else grid_dist_km(g, entry$cells)
  } else NULL
  lv <- if (model == "RR") NA_real_ else levels
  prob <- matrix(0, length(g$ei), length(t))
  status <- character(n_draws)
  metrics <- events <- vector("list", n_draws)
  t_start <- proc.time()[["elapsed"]]
  timed_out <- FALSE
  for (i in seq_len(n_draws)) {
    if (proc.time()[["elapsed"]] - t_start > time_limit) {
      timed_out <- TRUE
      status[i:n_draws] <- "timeout"
      break
    }
    pr <- draws[i, , drop = FALSE]
    cells_i <- if (model == "SLG") integer(0) else if (!is.null(entry)) entry$cells else rnd$entries[[i]]
    d_i <- if (model == "SLG") NULL else if (!is.null(fixed_d)) fixed_d else grid_dist_km(g, cells_i)
    res <- switch(model,
      SLG = {
        dens <- core_slg(g, pr$n0, pr$lmax, t)
        list(status = "ok", occ = lapply(levels, function(l) dens > l),
             years = sapply(levels, function(l) core_slg_years(g, pr$n0, pr$lmax, l)))
      },
      RR = {
        y <- core_radial_years(d_i, g$ei, pr$RR)
        list(status = "ok", occ = list(sapply(t, function(tt) y <= tt)), years = matrix(y, ncol = 1))
      },
      DK = {
        init <- ifelse(is.na(g$ei), NA, 0)
        init[cells_i] <- pr$init
        r <- do.call(core_dispk, c(list(g = g, init = init, lmax = pr$lmax, u = pr$u, p = pr$p, t = t, levels = levels),
                                   dk_args))
        if (!identical(r$status, "ok")) list(status = r$status) else
          list(status = "ok", occ = lapply(levels, function(l) r$density > l), years = r$years)
      })
    status[i] <- res$status
    if (res$status != "ok") next
    hi <- if (model == "RR") 1 else match(headline, levels)
    occ_h <- res$occ[[hi]]
    occ_h[is.na(occ_h)] <- FALSE
    prob <- prob + occ_h
    metrics[[i]] <- do.call(rbind, lapply(seq_along(lv), function(j) {
      occ <- res$occ[[j]]
      do.call(rbind, lapply(seq_along(t), function(k) {
        o <- occ[, k]
        o[is.na(o)] <- FALSE
        o <- o & suit
        data.frame(draw = i, t = t[k], level = lv[j], n_cells = sum(o), area_km2 = sum(g$area_km2[o]),
                   pct_suitable = 100 * sum(g$area_km2[o]) / sum(g$area_km2[suit]),
                   target_area_km2 = sum(g$area_km2[o & g$target]),
                   target_pct_suitable = if (any(suit & g$target)) 100 * sum(g$area_km2[o & g$target]) /
                     sum(g$area_km2[suit & g$target]) else NA_real_,
                   max_distance_km = if (!is.null(d_i) && any(o)) max(d_i[o]) else NA_real_)
      }))
    }))
    events[[i]] <- do.call(rbind, lapply(seq_along(lv), function(j) {
      y <- res$years[, j]
      tt <- tte_vec(y, g, pct, horizon)
      dd <- tte_vec(y, g_dom, pct, horizon)
      out <- data.frame(draw = i, level = lv[j], target_first_year = tt$first_year)
      for (x in pct) out[[paste0("target_year_pct_", x)]] <- tt[[paste0("year_pct_", x)]]
      for (x in pct) out[[paste0("domain_year_pct_", x)]] <- dd[[paste0("year_pct_", x)]]
      out
    }))
  }
  n_ok <- sum(status == "ok")
  add_id <- function(df) if (is.null(df) || !length(id)) df else cbind(as.data.frame(id, stringsAsFactors = FALSE),
                                                                       model = model, df)
  draws$draw <- seq_len(n_draws)
  draws$status <- status
  structure(list(draws = add_id(draws), metrics = add_id(do.call(rbind, metrics)), events = add_id(do.call(rbind, events)),
                 prob = if (n_ok) prob / n_ok else prob, n_ok = n_ok,
                 settings = list(model = model, t = t, levels = lv, headline = if (model == "RR") NA else headline,
                                 n_draws = n_draws, method = method, seed = seed, key = key, seed_used = s, pct = pct,
                                 priors = priors, id = id, horizon = horizon, timed_out = timed_out,
                                 elapsed_s = proc.time()[["elapsed"]] - t_start)),
            class = "sg_mc")
}

#' @export
print.sg_mc <- function(x, ...) {
  st <- x$settings
  cat(sprintf("<sg_mc> %s, %d draws (%s), %d ok, horizon %d yr, seed %d\n", st$model, st$n_draws, st$method, x$n_ok,
              st$horizon, st$seed_used))
  if (x$n_ok < st$n_draws) print(table(x$draws$status))
  invisible(x)
}

# Evaluate expr with a temporary seed; restores the caller's RNG state (or its absence)
with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had) assign(".Random.seed", old, envir = globalenv()) else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  expr
}
