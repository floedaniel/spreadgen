# spreadgen

Generic pest spread models of Robinet et al. (2012, *PLoS ONE* 7: e43366) for pest risk screening, ported to
modern R (`terra`) and run as Monte Carlo over PERT-distributed parameters:

| Model | Function | Parameters |
|---|---|---|
| SLG, simple logistic growth (paper: C) | `spread_slg()`, `core_slg()` | `lmax`, `n0` |
| RR, radial range expansion (B) | `spread_radial()`, `core_radial_years()` | `RR` (km/yr) |
| Rand-RR, random cells inside the radial range (A + B) | `spread_radialrand()` | `RR`, `n0`, `r` |
| DK, 2Dt dispersal kernel + logistic growth (D) | `spread_dispk()` (FFT), `spread_dispk_reference()` (literal port) | `u` (= RR), `p`, `lmax`, initial density |

LG-Econ (model A) is not implemented (no economic raster). Original code and inputs: DANS
doi:10.17026/DANS-2XV-Z3DX (CC0-1.0); `data-raw/download_original.R` fetches them, `data-raw/run_v21_baseline.R`
runs the original as the regression baseline. Deviations from the original: `inst/NOTES-deviations.md`;
kernel benchmarks: `inst/benchmarks-dispk.md`; worked example: `vignettes/spreadgen.Rmd`.

**Outputs are scenarios conditional on the assumptions (suitability, parameters, entry points), not predictions.**
The kernel shape `p` (nu) and `lmax` are weakly identifiable; GI from SDM suitability is an assumption (CLIMEX GI is
a population growth index); empirical front speeds include human-assisted spread.

## Main pieces

- `prep_suitability()` (SDM suitability → suitable mask at `threshold`, GI = 100 × suitability) or `prep_climex()`;
  `as_grid()` converts to a plain-data grid (vectors only) with a `target` mask (the reporting area).
- Density levels `c(0.1, 1, 10)` % of K (headline 1 %) for SLG and DK; first arrival from RR (`ceil(d / RR)`,
  `reach = "disc"` or `"land"`).
- `mc_spread()`: Latin hypercube / random draws, seeds per key (`seed_from_key()`), metrics per draw × t × level,
  time-to-event in the target (first cell, x % of its suitable area; "never within horizon" when censored),
  per-cell probability. `summarise_spread()` gives q05/q25/q50/q75/q95 and the mean.
- `run_batches()`: species-level parallelism on Windows (PSOCK workers via `future::multisession`, installed
  package, one thread per worker, FAILED rows, one sequential retry per crashed batch).

## Install

From the folder that contains `spreadgen/`:

    R CMD INSTALL spreadgen

Parallel workers load the installed package, so reinstall after every change (the batch driver checks that the
installed version matches and can reinstall itself).

## Tests

    Rscript -e "devtools::test('spreadgen')"
    set SPREADGEN_SLOW_TESTS=true && Rscript -e "devtools::test('spreadgen', filter = 'dispk')"   :: full 18-year V21 regression (~3 min)

## Long runs (the batch driver `4_scripts/spread_robinet.R`)

- Run with `Rscript` (or `4_scripts/run_spread_robinet.bat`), not inside the RStudio session: a long session can
  be interrupted by the IDE, and PSOCK workers are separate R processes anyway.
- Disable sleep / hibernate for the duration (Windows: Settings > System > Power), or run
  `powercfg /change standby-timeout-ac 0` before and restore afterwards.
- Optionally pause OneDrive sync during the run: per-species outputs are written to a local folder
  (`LOCAL_DIR`, default `C:/tmp/spreadgen`) and copied into the species folders once per batch.
- Every setting is in chapter 0 of `spread_robinet.R`; the log goes to `4_scripts/logs/`.
