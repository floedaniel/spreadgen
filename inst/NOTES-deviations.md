# Deviations from the original code (Prog_SpreadModule_V21-web.R)

Every intended difference between `spreadgen` and the original V21 script is listed here, with its reason and its
measured effect on the D. virgifera regression (fixture `tests/testthat/fixtures/v21_dvv_baseline.rds`).

## Observations on the original (M0)

- **O1. Cell area.** Area is a flat-earth approximation: `(40000/360)·dx·cos(lat) × (40000/360)·dy` km² with
  dx, dy in degrees (Earth circumference 40,000 km). It is applied as the landing-cell area in `dispk()`.
- **O2. Distances.** Great circle on a sphere of radius 20000/π ≈ 6366.2 km; the cosine argument is rounded to 7
  digits before `acos` (≈ 3 km resolution at very short distance; irrelevant at 0.5°).
- **O3. Sphere correction.** `dispk()` multiplies the 2Dt density by `1/s`, `s = sin(R/ρ)/(R/ρ)`, ρ = 20000/π.
- **O4. Sources.** `dispk()` loops over every suitable cell (EI > 0) or presence cell, whether or not it is occupied;
  unoccupied sources add 0. The reset `popdisp[EI == 0] <- 0` runs inside the source loop; it is equivalent to one
  reset after the loop.
- **O5. Two grids.** SLG / radial use the CLIMEX point list (`tab`); `dispk()` uses the full rectangular lon/lat grid
  (`vectxy`, sea cells EI = 0). Counts `ntot` differ between the two.
- **O6. Snapping.** Entry points (radial) and presence points (`dispk`) are snapped to the nearest cell by squared
  degree differences (not great-circle distance); radial snaps to suitable cells only, `dispk` to any grid cell.
- **O7. Invaded threshold.** `dispk()` counts a cell as invaded when density > 1 % of K (strict).

- **O8. Baseline run (2026-10-19, R 4.5.3, raster 3.x; `data-raw/run_v21_baseline.R`).** Reproduces the published
  D. virgifera areas exactly: SLG t = 8 (> 1 % of K) 59,790 km²; radial t = 18 3,578,880 km² (56.3 % of 3104
  suitable cells, "more than half"); DK t = 18 2,010,382 km² (976 cells, 31.4 %). DK is 0 cells above 1 % until
  t = 16 (densities 1e-6..0.5 %), then 90 / 399 / 976 cells at t = 16 / 17 / 18: invaded-area metrics for DK are
  threshold- and horizon-sensitive (cf. paper Tables 2-3). Runtime: DK 590 s for 18 years (≈ 33 s/yr).
- **O9. Unresolved vs documentation.** The tutorial quotes 17,298 DK grid cells and 3,326 suitable DK cells for the DD
  version; this run has 20,995 and 3,104 (the tutorial predates V21 / used another input extent). Population totals
  (`sum(N/100) × 6.3e10`) give 1.0e11 (SLG, t = 8) and 3.1e12 (DK, t = 18) vs 3.2e10 and 1.6e12 in the paper: the
  paper's totals likely use another aggregation (e.g. host-area weighting); areas, which the regression uses, match.

## Deviations (spreadgen)

| # | Deviation | Reason | Effect on DvV regression |
|---|---|---|---|
| D1 | GI = suitability × 100 without rescaling to the PRA maximum (`gi_rescale = "none"`, default); `"domain"` reproduces the original rescaling | User decision 2026-10-19: suitability is a probability, rescaling to the domain maximum would change the growth rate with the extent | n/a (regression uses CLIMEX GI, rescaled as the original) |
| D2 | LG-Econ (model A) not implemented | User decision 2026-10-19: no economic raster | n/a |
| D3 | GI set to 0 outside suitable cells (original masks GI by host / elevation only) | Clarity; every model zeroes density where EI = 0 before or after growth | None: SLG identical to 1e-12 over t = 1..18 |
| D4 | Cell area default `terra::cellSize()` (WGS84 ellipsoid, or planar on projected grids); `area_method = "v21"` gives the original sphere approximation (O1) | Exact areas; equal-area grids in production | Ellipsoid vs v21 up to 1.02 % per cell on the 0.5° grid; regression tests use `"v21"` (published SLG area 59,790 km² reproduced) |
| D5 | Entry snapping by great-circle (lon/lat, WGS84) or Euclidean (projected) distance instead of squared degree differences (O6) | Correct nearest cell at high latitudes | Belgrade snaps to the same cell (20.25 E, 44.75 N) |
| D6 | Two density thresholds, `reach` (default 1e-6 % of K) and `abundant` (default 1 %, the original's invaded threshold, strict `>`) | User decision 2026-10-19: DK / SLG invaded counts are threshold- and horizon-sensitive (O8) | `abundant` = original `ninv` |
| D7 | Logistic evaluated as `N0 / (N0/100 + (1 - N0/100) e^(-rt))` | Same value; no overflow (NaN) for large `r t` | Identical to 1e-12 |
| D8 | Population totals not reported | User decision 2026-10-19 (O9); host-weighted K remains possible via `k_carrying(area, host_raster, density)` | n/a |
| D9 | Original SLG class counts (`sum0..sum100`, with the original's overlapping boundary at 25) only in `slg_classes_v21()` | Regression; the summary uses the two thresholds | Identical counts t = 1..18 |
| D10 | RR distance modes: `"sphere"` (default; the original great circle, radius 20000/π km, cosine clamped to [-1, 1], no rounding), `"v21"` (cosine rounded to 7 digits before `acos`, inclusive `R <= RR*t`), `"grid"` (`terra::distance`: WGS84 geodesic / Euclidean on projected grids, fastest with many entry cells) | The 7-digit rounding shifts short distances by up to ~2 km (1° → −1.8 km); `"v21"` kept for exact regression | Jaccard vs V21, t = 1..18: `v21` 1, `sphere` 1, `grid` ≥ 0.986 (≤ 6 cells) |
| D11 | `reach = "land"`: distance as the shortest path through cells with data (`terra::gridDist`), sea / no-data cells are barriers. **Not in Robinet et al. (2012)**, which uses straight-line discs (`reach = "disc"`, default) | Comparability with the route-over-land screen of the sibling scripts | n/a (option; not used in the regression) |
| D12 | RR also returns `years_to_reach = ceil(d / RR)` (d = distance to the nearest entry cell; entry cells 0) | Time-to-event outputs (user decision 2026-10-19); equals the first integer `t` with `d <= RR*t` | Consistent with the invaded stacks (tested) |
| D13 | Rand-RR draws one random order of cells per run, so the invaded set is nested in time; it uses the caller's RNG stream (no seed argument; `mc_spread()` seeds) | Original drew independently per call without seed control (issue 11) | Counts follow the original formula (tested); cell sets are random |

The disc is the union of discs around the entry cells only (as V21: newly invaded cells do not start discs).
Entry snapping ties (V21: first cell by squared degree distance) are not tested.

(Extended as milestones add code.)
