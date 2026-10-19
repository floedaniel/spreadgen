# Original generic spread module (Robinet et al. 2012)

Source: Robinet C, Kehlenbeck H, Kriticos DJ, et al. (2012) A suite of models to support the quantitative assessment
of spread in pest risk analysis. *PLoS ONE* 7(10): e43366. doi:10.1371/journal.pone.0043366

Code and inputs: DANS Data Station Life Sciences, doi:10.17026/DANS-2XV-Z3DX. Licence: **CC0-1.0** (public domain).

The files are not stored in git. Re-create them from the package root with

    Rscript data-raw/download_original.R        # PDFs, program versions, D. virgifera inputs (~10 MB)
    Rscript data-raw/download_original.R all    # every case study (incl. a 188 MB ESRI grid); no installers / old zips

The file list with DANS ids and folders is `data-raw/dans_manifest.csv` (file names repeat across case-study folders;
the `dir` column is DANS `directoryLabel`).

*Diabrotica virgifera virgifera* (case study 1) = folder `Code for case studies/1-Diabrotica virgifera virgifera`:
`ClimexOutput.txt` (23576), `habitat.tif` (23513), `econ.tif` (23562), `presence.txt` (23460), `europeL.txt` (23538),
`Prog_SpreadModule_V21-web.R` (23542), `useful commands.r` (23557).

The patched baseline used for regression tests is `data-raw/run_v21_baseline.R`; it does not modify these files.
