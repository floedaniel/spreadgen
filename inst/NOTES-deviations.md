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

(Extended as milestones add code.)
