# Parameter priors from Robinet et al. (2012), Table 1: species values and taxon-level defaults (low confidence).

#' Robinet et al. (2012) Table 1 parameter values
#'
#' Ranges as published; `NA` = not given ("n/a"). `Pmax` = carrying capacity per cell (units differ per species),
#' `lmax` = maximum yearly multiplication factor, `RR` = radial rate c (km/yr), `u` = kernel scale (km),
#' `p` = kernel shape (the paper's nu).
#' @return Data frame, one row per species.
#' @export
robinet_table1 <- function() {
  data.frame(
    species = c("Diabrotica virgifera virgifera", "Anoplophora chinensis", "Anoplophora glabripennis", "Eichhornia crassipes",
                "Meloidogyne enterolobii", "Bursaphelenchus xylophilus / Monochamus", "Gibberella circinata"),
    group = c("insect", "insect", "insect", "plant", "nematode", "nematode", "pathogen"),
    r_min = c(0.33, NA, NA, NA, NA, 0.27, NA), r_max = c(0.33, NA, NA, NA, NA, 0.27, NA),
    RR_min = c(80, 1, 1.5, 30, 10, 6, NA), RR_max = c(80, 2, 3, 100, 30, 35, NA),
    Pmax_min = c(6.3e10, 1.9e7, 7.58e6, 7e7, 6.2e14, 14.69, 5.21e8), Pmax_max = c(6.3e10, 1.9e8, 1.5e7, 7e7, 6.2e14, 14.69, 5.21e8),
    lmax_min = c(40, 6, 5, 30, 9.7, 8.76, 2.72), lmax_max = c(40, 6, 5, 30, 9.7, 8.76, 2.72),
    u_min = c(80, 2, 1.5, 30, 10, 35, 1), u_max = c(80, 2, 3, 70, 30, 35, 1),
    p_min = c(5, 10, 30, 10, 10, 5, 2), p_max = c(5, 50, 50, 50, 50, 5, 2),
    stringsAsFactors = FALSE)
}

#' Taxon group for the Table 1 defaults
#'
#' @param kingdom,phylum Taxonomy (e.g. from the master species list).
#' @return `"insect"` (Arthropoda), `"nematode"` (Nematoda), `"plant"` (Plantae) or `"pathogen"` (fungi, oomycetes,
#'   bacteria, viruses and everything else).
#' @export
taxon_group <- function(kingdom = NA, phylum = NA) {
  k <- tolower(as.character(kingdom))
  ph <- tolower(as.character(phylum))
  ifelse(!is.na(ph) & ph == "arthropoda", "insect",
         ifelse(!is.na(ph) & ph == "nematoda", "nematode",
                ifelse(!is.na(k) & k == "plantae", "plant", "pathogen")))
}

#' Taxon-level default priors from Robinet et al. (2012) Table 1 (low confidence)
#'
#' Pools the Table 1 species of the group: PERT(min, mode, max) on the log scale with min / max = the smallest /
#' largest published value and mode = the geometric mean of the species' values; a group with one value gets +-25 %.
#' Used for `lmax` and `p` (no other source), and for `RR` only as the last fallback.
#' @param group One of `taxon_group()`'s values.
#' @param table Table 1 ([robinet_table1()]).
#' @return Named list of `sg_prior` (`lmax`, `p`, and `RR` when the group has values); `source` and
#'   `confidence = "low"` set.
#' @export
taxon_priors <- function(group, table = robinet_table1()) {
  rows <- table[table$group == group, , drop = FALSE]
  if (!nrow(rows)) stop("no Table 1 species in group '", group, "'")
  src <- paste0("Robinet et al. 2012 Table 1, ", group, " default (", paste(rows$species, collapse = "; "), ")")
  pool <- function(lo, hi, lower = 0) {
    v <- c(lo, hi)
    v <- v[!is.na(v)]
    if (!length(v)) return(NULL)
    mode <- exp(mean(log(sqrt(lo * hi)), na.rm = TRUE))
    mn <- min(v)
    mx <- max(v)
    if (mx / mn < 1 + 1e-9) { mn <- mode * 0.75; mx <- mode * 1.25 }
    prior_pert(max(mn, lower), min(max(mode, max(mn, lower)), mx), mx, log = TRUE, source = src, confidence = "low")
  }
  out <- list(lmax = pool(rows$lmax_min, rows$lmax_max, lower = 1.0001), p = pool(rows$p_min, rows$p_max, lower = 1.05),
              RR = pool(rows$RR_min, rows$RR_max))
  out[!vapply(out, is.null, logical(1))]
}
