# Inputs for the batch driver (not part of the R package)

`spread_robinet_priors.xlsx` is the single input file with all parameter values (min / most likely / max, PERT)
used by the driver `spread_robinet.R` (VKM plant health data project, `4_scripts/`). One row per species x
parameter; the sheet `README` inside the file explains the columns.

- Built once by the driver from `spread_plot.R` front speeds (RR "empirical total"), the trait table
  `spread_trait_groups.xlsx` (RR "natural"), Robinet et al. (2012) Table 1 taxon defaults (`lmax`, `p`) and default
  starting densities. Afterwards this file is the only source: edit the values here.
- Species not yet in the file are appended on the next run; existing rows are never changed.
- Changing a species' rows makes `RESUME` rerun that species.

This folder is excluded from package builds (`.Rbuildignore`).
