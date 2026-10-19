# DK benchmarks (M3, 2026-10-19)

Script: `data-raw/benchmark_dispk.R`; raw results `data-raw/benchmark_dispk_results.csv` (A) and
`data-raw/benchmark_dispk_cutoff.csv` (B). Machine: the project laptop (16 physical / 22 logical cores, 63 GB), one
core, R 4.5.3, base `fft()`. Grid: land in Europe north of 45 N, EPSG:3035 (same construction as the sibling scripts),
synthetic suitability 0.5 on land; 5 entry cells near Hamburg at 10 % of K; 20 years per draw.

## Regression

| Check | Result |
|---|---|
| Reference implementation vs V21 (0.5° lon/lat, v21 areas, t = 1..18) | max abs 1.1e-13, rel 1.2e-14; 976 cells, 2,010,382 km² at t = 18; 196 s (V21: 590 s) |
| FFT vs reference, same projected grid | ≤ 1e-14 rel (unit test, small grid); 8.6e-12 rel on the 25 km Europe grid |
| FFT (25 km LAEA) vs V21 (0.5° lon/lat), D. virgifera t = 18 | 978,773 vs 2,010,382 km² (−51 %); with the same founder number per entry cell (N0 rescaled to the 625 km² cell's K) 1,747,542 km² (−13 %). N0 is % of the cell's K, so a smaller entry cell holds fewer founders; DK invaded area at t = 16..18 is in its exponential phase (V21: 194k / 840k / 2,010k km²), so a fraction of a year shifts it strongly |

## A. Time per draw (20 years, one core)

| Grid | Cells | DvV-like u 80, p 5 (auto = full grid) | fat tail u 20, p 2 (auto) | local u 15, p 50 (auto R_c 155 km) | DvV-like, k = 10 crop |
|---|---|---|---|---|---|
| 20 km | 30,653 | 0.8 s | skipped (u < 2 cells) | skipped | 0.3 s |
| 10 km | 121,200 | 2.6 s | 2.4 s | skipped | 0.3 s |
| 5 km | 482,586 | 11.5 s | 15.0 s | 4.6 s | 1.2 s |

Peak memory ≤ 0.42 GB. Thin-tailed kernels get a small automatic cutoff; for p ≤ ~5 the conservative automatic
cutoff is the whole grid (FFT 1200 × 1620 at 5 km).

## B. Cutoff sensitivity (10 km; Norway time-to-event, years)

| Kernel | k | R_c km | tail mass | first reach | 10 % reached | 50 % reached | first abundant |
|---|---|---|---|---|---|---|---|
| DvV-like | 3 / 5 / 10 / 20 / auto | 240 / 400 / 800 / 1600 / grid | 0.13 / 0.028 / 0.0023 / 1.5e-4 / 4e-8 | 3 / 2 / 1 / 1 / 1 | 3 / 2 / 1 / 1 / 1 | 6 / 3 / 2 / 1 / 1 | 6 / 5 / 5 / 5 / 5 |
| fat tail | 3 / 5 / 10 / 20 / auto | 60 / 100 / 200 / 400 / grid | 0.43 / 0.27 / 0.14 / 0.07 / 9e-5 | never / 11 / 3 / 2 / 1 | never / 16 / 6 / 2 / 1 | never / never / 11 / 5 / 1 | never / never / 13 / 9 / 8 |

Truncation changes time-to-reach by years for fat tails (and even `k = 10` is not enough); the "abundant" (1 % of K)
times are less sensitive. Default `cutoff = "auto"` (conservative, summed over all suitable cells).

## Recommendation for Monte Carlo

- Kernel model on the 5 km grid: ~5-15 s per draw → 200 Latin-hypercube draws ≈ 15-50 min on one core, ~3-6 min on
  10 workers (parallel over draws). 10 km (2.5 s/draw) for quick screening.
- 200 LHS draws for DK (u, p, lmax, initial density); 1000 plain draws for SLG and RR (closed form, < 0.1 s each).
- u < 2 cells (10 km at 5 km resolution) skips DK (`spread below grid`), as decided.
- The reach threshold (1e-6 % of K) is passed at long range in year 1 by the deterministic kernel tail (Norway
  reached in year 1 from Hamburg for u = 80 or p = 2): report "abundant" alongside "reach".
