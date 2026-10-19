# ---- priors and sampling ----

test_that("PERT: bounds, mode, mean; degenerate prior", {
  set.seed(1)
  x <- rpert(2e5, 2, 5, 14)
  expect_true(all(x >= 2 & x <= 14))
  expect_equal(mean(x), (2 + 4 * 5 + 14) / 6, tolerance = 5e-3)
  d <- density(x)
  expect_equal(d$x[which.max(d$y)], 5, tolerance = 0.1)
  expect_equal(qpert(c(0.1, 0.9), 3, 3, 3), c(3, 3))
  expect_error(prior_pert(3, 2, 4))
  expect_error(prior_pert(0, 1, 2, log = TRUE))
})

test_that("Latin hypercube: one draw per stratum in every dimension", {
  set.seed(2)
  U <- lhs_unif(50, 3)
  for (j in 1:3) expect_equal(sort(floor(U[, j] * 50)), 0:49)
})

test_that("draw_params respects bounds, log scale and value / range helpers", {
  set.seed(3)
  pr <- list(RR = prior_pert(1, 10, 100, log = TRUE), p = prior_pert(2, 5, 50))
  d <- draw_params(pr, 400)
  expect_true(all(d$RR >= 1 & d$RR <= 100 & d$p >= 2 & d$p <= 50))
  expect_equal(median(log(d$RR)), log(10), tolerance = 0.1)
  v <- prior_from_value(80)
  expect_equal(c(v$min, v$mode, v$max), c(60, 80, 100))
  r <- prior_from_range(1, 4)
  expect_equal(r$mode, 2)
})

test_that("seed_from_key is deterministic and key-specific", {
  expect_equal(seed_from_key(1, "DIABVI"), seed_from_key(1, "DIABVI"))
  expect_false(seed_from_key(1, "DIABVI") == seed_from_key(1, "PAYSAR"))
  expect_false(seed_from_key(1, "DIABVI") == seed_from_key(2, "DIABVI"))
  expect_equal(seed_from_key(7), 7L)
})

test_that("Table 1 taxon priors: low confidence, bounds from the published values", {
  tp <- taxon_priors("insect")
  expect_equal(c(tp$lmax$min, tp$lmax$max), c(5, 40))
  expect_equal(c(tp$p$min, tp$p$max), c(5, 50))
  expect_equal(tp$lmax$confidence, "low")
  pa <- taxon_priors("pathogen")
  expect_equal(c(pa$lmax$min, pa$lmax$max), 2.72 * c(0.75, 1.25))
  expect_null(pa$RR)
  expect_equal(taxon_group(c("Animalia", "Animalia", "Plantae", "Fungi"), c("Arthropoda", "Nematoda", NA, "Ascomycota")),
               c("insect", "nematode", "plant", "pathogen"))
})

test_that("read_empirical_priors: usable rows, floor flag, reasons, rr_source", {
  x <- data.frame(eppo = c("A", "B", "C", "D", "E", "F"), front_km_yr = c(10, NA, 5, 8, 8, -1),
                  front_lo = c(6, NA, 2, 0, -3, -5), front_hi = c(15, NA, 9, 12, 12, 4),
                  flags = c("", "", "short window (6 yr)", NA, "", ""))
  e <- read_empirical_priors(x)
  expect_equal(e$usable, c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE))
  expect_equal(e$floor_needed, c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE))
  expect_match(e$reason[3], "short window")
  expect_equal(unique(e$rr_source), "empirical front speed")
  u <- read_empirical_priors(x, nonpositive_lo = "unusable")
  expect_equal(u$usable, c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
})

test_that("rr_floor: trait-table lower bound, else 0.1 x median; positive bounds kept", {
  f <- rr_floor(c(-3, 0, 2, -1), c(8, 8, 5, 10), c(1.5, NA, 1, 20))
  expect_equal(f$rr_min, c(1.5, 0.8, 2, 1))
  expect_equal(f$floor_applied, c(TRUE, TRUE, FALSE, TRUE))
  expect_equal(f$floor_source, c("trait table", "0.1 x median", NA, "0.1 x median"))
})

# ---- mc_spread ----

mc_grid <- function() {
  s <- synth_suit(20, 30) * 0.6 + 0.4
  s[3:6, 20:23] <- 0.05
  p <- prep_suitability(s, 0.1)
  tgt <- s * 0
  tgt[11:20, 16:30] <- 1                                    # target: south-east quarter
  tgt <- terra::ifel(tgt == 1, 1, NA)
  as_grid(p, target = tgt)
}

test_that("mc_spread is reproducible, seeded by key, and restores the caller's RNG", {
  g <- mc_grid()
  pr <- list(RR = prior_pert(1, 2, 4, log = TRUE))
  set.seed(99)
  before <- .Random.seed
  a <- mc_spread(g, "RR", pr, 12, t = c(2, 5, 10), entry = list(cells = 35L), seed = 1, key = "X")
  expect_identical(.Random.seed, before)
  b <- mc_spread(g, "RR", pr, 12, t = c(2, 5, 10), entry = list(cells = 35L), seed = 1, key = "X")
  expect_identical(a$draws, b$draws)
  expect_identical(a$metrics, b$metrics)
  c2 <- mc_spread(g, "RR", pr, 12, t = c(2, 5, 10), entry = list(cells = 35L), seed = 1, key = "Y")
  expect_false(identical(a$draws$RR, c2$draws$RR))
})

test_that("mc_spread RR: metrics agree with the core; probability in [0, 1]; target events", {
  g <- mc_grid()
  m <- mc_spread(g, "RR", list(RR = prior_pert(1, 2, 4, log = TRUE)), 10, t = c(3, 8), entry = list(cells = 35L),
                 id = list(species = "TEST", scenario = "natural"))
  expect_equal(m$n_ok, 10)
  d <- grid_dist_km(g, 35L)
  y <- core_radial_years(d, g$ei, m$draws$RR[4])
  r4 <- m$metrics[m$metrics$draw == 4 & m$metrics$t == 8, ]
  expect_equal(r4$n_cells, sum(y <= 8, na.rm = TRUE))
  expect_true(all(m$prob >= 0 & m$prob <= 1))
  ev <- m$events[m$events$draw == 4, ]
  suit_t <- g$target & !is.na(g$ei) & g$ei == 1
  first <- min(y[suit_t])
  expect_equal(ev$target_first_year, if (first <= 8) first else Inf)   # censored at the horizon (t = 8)
  expect_equal(unique(m$metrics$species), "TEST")
  expect_true(all(is.na(m$metrics$level)))                  # RR: arrival, no density level
})

test_that("mc_spread SLG and DK: levels, skip below grid, sampled entries respect weights", {
  g <- mc_grid()
  s <- mc_spread(g, "SLG", list(lmax = prior_pert(3, 5, 10, log = TRUE), n0 = prior_pert(1e-5, 1e-4, 1e-3, log = TRUE)),
                 8, t = c(5, 10))
  expect_equal(sort(unique(s$metrics$level)), c(0.1, 1, 10))
  expect_equal(nrow(s$events), 8 * 3)
  dk_pr <- list(RR = prior_pert(1.5, 3, 4, log = TRUE), p = prior_pert(3, 5, 10), lmax = prior_pert(3, 5, 8, log = TRUE),
                init = prior_pert(1, 10, 50, log = TRUE))
  k <- mc_spread(g, "DK", dk_pr, 6, t = c(2, 4), entry = list(cells = c(35L, 36L)), dk_args = list(cutoff = "none"))
  expect_true(all(k$draws$status[k$draws$u < 2] == "spread below grid"))
  expect_true(all(k$draws$status[k$draws$u >= 2] == "ok"))
  expect_equal(k$draws$u, k$draws$RR)                       # u = RR rule
  w <- rep(0, length(g$ei))
  w[c(500, 510)] <- c(1, 3)
  r <- mc_spread(g, "RR", list(RR = prior_pert(1, 2, 3)), 20, t = 2, entry_sampler = list(weights = w, n = 1))
  first <- r$metrics$n_cells
  expect_true(all(first >= 1))
})

test_that("summarise_spread and fmt_interval", {
  g <- mc_grid()
  m <- mc_spread(g, "RR", list(RR = prior_pert(1, 2, 4, log = TRUE)), 20, t = c(3, 8), entry = list(cells = 35L),
                 id = list(species = "TEST"))
  s <- summarise_spread(m)
  expect_true(all(c("species", "model", "t", "metric", "mean", "q05", "q25", "q50", "q75", "q95") %in% names(s$metrics)))
  row <- s$metrics[s$metrics$t == 8 & s$metrics$metric == "n_cells", ]
  expect_equal(row$q50, median(m$metrics$n_cells[m$metrics$t == 8]))
  expect_true("share_never" %in% names(s$events))
  expect_equal(fmt_interval(c(1, 2, 3, 4, 100)), "3 (1-90)")
  expect_equal(fmt_interval(c(Inf, Inf, 2)), "never within horizon")
  expect_match(fmt_interval(c(2, 3, Inf, Inf, Inf, 4, 5, 6, 7, 8)), "never within horizon\\)$")
})

test_that("plots build", {
  g <- mc_grid()
  m <- mc_spread(g, "RR", list(RR = prior_pert(1, 2, 4, log = TRUE)), 5, t = c(3, 8), entry = list(cells = 35L),
                 id = list(label = "natural"))
  expect_s3_class(plot_probability(g, m$prob[, 2], "p"), "ggplot")
  expect_s3_class(plot_events(m$events, horizon = 8), "ggplot")
  expect_s3_class(plot_metric(m$metrics), "ggplot")
})

test_that("sg_workers respects caps", {
  expect_equal(sg_workers(max_workers = 3, avail_gb = 100), min(3L, max(1L, parallel::detectCores(logical = FALSE) - 2L)))
  expect_equal(sg_workers(avail_gb = 1, mem_per_worker_gb = 0.5), min(2L, max(1L, parallel::detectCores(logical = FALSE) - 2L)))
  expect_equal(sg_workers(avail_gb = 0.1), 1L)
})

test_that("mc_spread time_limit stops between draws and flags the run", {
  g <- mc_grid()
  pr <- list(RR = prior_pert(1, 2, 4, log = TRUE))
  a <- mc_spread(g, "RR", pr, 6, t = c(2, 5), entry = list(cells = 35L), seed = 1, key = "X", time_limit = -1)
  expect_true(a$settings$timed_out)
  expect_equal(a$n_ok, 0)
  expect_true(all(a$draws$status == "timeout"))
  b <- mc_spread(g, "RR", pr, 6, t = c(2, 5), entry = list(cells = 35L), seed = 1, key = "X")
  expect_false(b$settings$timed_out)
  expect_equal(b$n_ok, 6)
})
