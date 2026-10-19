# Download the original Robinet et al. (2012) spread-module code and inputs from DANS (CC0-1.0).
# doi:10.17026/DANS-2XV-Z3DX. Large files (ESRI grids, R installers, old package zips) are not fetched by default.
# Files land in inst/original/<directoryLabel>/<filename>; the default set is stored in git, further files are ignored.
#
# Usage (from the package root): Rscript data-raw/download_original.R [all]

DANS_API <- "https://lifesciences.datastations.nl/api/access/datafile/"
OUT      <- file.path("inst", "original")
MAX_MB   <- 20                     # skip files above this size unless 'all'

args <- commandArgs(trailingOnly = TRUE)
get_all <- "all" %in% args

man <- utils::read.csv(file.path("data-raw", "dans_manifest.csv"), stringsAsFactors = FALSE)
man$dir[is.na(man$dir)] <- ""

# Default set: top-level PDFs, the three program versions, D. virgifera case-study inputs
dvv <- "Code for case studies/1-Diabrotica virgifera virgifera"
keep <- (man$dir == "" & grepl("\\.pdf$", man$name)) |
  man$name %in% c("Prog_SpreadModule_V21-web.R", "Prog_SpreadModule_metric_V08-web.R", "Prog_SpreadModule_V21world-web.R") &
    man$dir %in% c(dvv, "Spread module code/metric version-Europe", "Spread module code/DD version-world") |
  man$dir == dvv
if (get_all) keep <- !grepl("\\.(exe|zip)$", man$name)
sel <- man[keep & (get_all | man$size < MAX_MB * 1e6), ]

for (i in seq_len(nrow(sel))) {
  dest <- file.path(OUT, sel$dir[i], sel$name[i])
  if (file.exists(dest) && file.size(dest) == sel$size[i]) next
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  utils::download.file(paste0(DANS_API, sel$id[i]), dest, mode = "wb", quiet = TRUE)
  ok <- file.size(dest) == sel$size[i]
  message(sprintf("%-6s %s %s", sel$id[i], if (ok) "ok  " else "SIZE MISMATCH", file.path(sel$dir[i], sel$name[i])))
}
