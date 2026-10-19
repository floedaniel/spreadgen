# spreadgen

Generic pest spread models of Robinet et al. (2012) for pest risk screening, ported to modern R (`terra`) and run
as Monte Carlo over PERT-distributed parameters. Introduction, models, main functions and a worked example:
`vignette("spreadgen")`.

> Robinet C, Kehlenbeck H, Kriticos DJ, Baker RH, Battisti A, Brunel S, Dupin M, Eyre D, Faccoli M, Ilieva Z,
> Kenis M, Knight J, Reynaud P, Yart A, van der Werf W (2012). A suite of models to support the quantitative
> assessment of spread in pest risk analysis. *PLoS ONE* 7(10): e43366. doi:10.1371/journal.pone.0043366

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

From this repository (e.g. `C:/Dev/spreadgen`):

    R CMD INSTALL .

or from GitHub: `remotes::install_github("floedaniel/spreadgen")`.

Parallel workers load the installed package, so reinstall after every change (the batch driver checks that the
installed version matches and can reinstall itself).

## Tests

    Rscript -e "devtools::test()"
    set SPREADGEN_SLOW_TESTS=true && Rscript -e "devtools::test(filter = 'dispk')"   :: full 18-year V21 regression (~3 min)

## Long runs (a batch driver, e.g. `spread_robinet.R` in the VKM plant health data project)

The package has no project paths. A driver script keeps the inputs and settings and points to this repository as
the package source. Notes for long runs on Windows:

- Run with `Rscript` (or a `.bat` file), not inside the RStudio session: a long session can
  be interrupted by the IDE, and PSOCK workers are separate R processes anyway.
- Disable sleep / hibernate for the duration (Windows: Settings > System > Power), or run
  `powercfg /change standby-timeout-ac 0` before and restore afterwards.
- Optionally pause OneDrive sync during the run: per-species outputs are written to a local folder
  (`LOCAL_DIR`, default `C:/tmp/spreadgen`) and copied into the species folders once per batch.
- In `spread_robinet.R` every setting is in chapter 0 (`SPREADGEN_SRC` = this repository); the log goes to
  `4_scripts/logs/`.
