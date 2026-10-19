# Original generic spread module (Robinet et al. 2012)

Article: Robinet C, Kehlenbeck H, Kriticos DJ, Baker RH, Battisti A, Brunel S, Dupin M, Eyre D, Faccoli M, Ilieva Z,
Kenis M, Knight J, Reynaud P, Yart A, van der Werf W (2012). A suite of models to support the quantitative assessment
of spread in pest risk analysis. *PLoS ONE* 7(10): e43366. doi:10.1371/journal.pone.0043366. PMID: 23056174;
PMCID: PMC3467266. The article (open access, CC BY) is `Robinet_2012_PLoSONE.pdf`.

Code, tutorial and case-study inputs: DANS Data Station Life Sciences, doi:10.17026/DANS-2XV-Z3DX, version 1.0.
Licence: **CC0-1.0** (public domain).

## Files in this folder

These files are stored in git (they are excluded from package builds by `.Rbuildignore`, so the installed package
does not carry them). They are used by `data-raw/run_v21_baseline.R` (regression baseline, which does not modify
them) and for `inst/extdata`.

| File | kB | DANS file id |
|---|---|---|
| `Case studies.pdf` | 5,126 | [23519](https://lifesciences.datastations.nl/file.xhtml?fileId=23519&version=1.0) |
| `Code for case studies/1-Diabrotica virgifera virgifera/ClimexOutput.txt` |   133 | [23576](https://lifesciences.datastations.nl/file.xhtml?fileId=23576&version=1.0) |
| `Code for case studies/1-Diabrotica virgifera virgifera/econ.tif` |   325 | [23562](https://lifesciences.datastations.nl/file.xhtml?fileId=23562&version=1.0) |
| `Code for case studies/1-Diabrotica virgifera virgifera/europeL.txt` |   482 | [23538](https://lifesciences.datastations.nl/file.xhtml?fileId=23538&version=1.0) |
| `Code for case studies/1-Diabrotica virgifera virgifera/habitat.tif` |   166 | [23513](https://lifesciences.datastations.nl/file.xhtml?fileId=23513&version=1.0) |
| `Code for case studies/1-Diabrotica virgifera virgifera/presence.txt` |     0 | [23460](https://lifesciences.datastations.nl/file.xhtml?fileId=23460&version=1.0) |
| `Code for case studies/1-Diabrotica virgifera virgifera/Prog_SpreadModule_V21-web.R` |    46 | [23542](https://lifesciences.datastations.nl/file.xhtml?fileId=23542&version=1.0) |
| `Code for case studies/1-Diabrotica virgifera virgifera/SpreadModule.RData` |     0 | [23453](https://lifesciences.datastations.nl/file.xhtml?fileId=23453&version=1.0) |
| `Code for case studies/1-Diabrotica virgifera virgifera/useful commands.r` |     1 | [23557](https://lifesciences.datastations.nl/file.xhtml?fileId=23557&version=1.0) |
| `Spread module code/DD version-world/Prog_SpreadModule_V21world-web.R` |    42 | [23504](https://lifesciences.datastations.nl/file.xhtml?fileId=23504&version=1.0) |
| `Spread module code/metric version-Europe/Prog_SpreadModule_metric_V08-web.R` |    33 | [23570](https://lifesciences.datastations.nl/file.xhtml?fileId=23570&version=1.0) |
| `Tutorial on the generic spread models.pdf` | 1,593 | [23540](https://lifesciences.datastations.nl/file.xhtml?fileId=23540&version=1.0) |

*Diabrotica virgifera virgifera* is case study 1 of the paper.

## The rest of the DANS dataset

The full dataset is 316 MB in 124 files (all case studies, ESRI grids, R installers, old package zips); one ESRI
grid file is 189 MB, above the GitHub file limit, so the rest is not stored here. The file list with DANS ids and
folders is `data-raw/dans_manifest.csv` (file names repeat across case-study folders; `dir` = DANS `directoryLabel`).
Download from the package root with

    Rscript data-raw/download_original.R        # the files above
    Rscript data-raw/download_original.R all    # every case study (no installers / old zips)

Files downloaded beyond the list above are git-ignored.
