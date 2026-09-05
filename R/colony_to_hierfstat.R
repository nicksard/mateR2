#' @title Convert Simulated COLONY Genotypes to hierfstat Format
#'
#' @description Reshapes the genotypes read by \code{\link{read_colony_dat}}
#'   into the data frame \pkg{hierfstat} expects: a grouping column followed by
#'   one column per locus, each genotype encoded as its two alleles concatenated
#'   into a single number. Because a mateR2 cohort comes with a known pedigree,
#'   the grouping can be the true full-sibling family, the paternal or maternal
#'   half-sibling family, or the cohort as a whole.
#'
#' @details
#' \strong{Why the grouping matters.} With one simulated cohort drawn from a
#' single allele-frequency pool there is no population structure to find, so
#' \code{level = "cohort"} is only useful for \eqn{H_o}, \eqn{H_s} and departures
#' from Hardy-Weinberg. The informative use is \code{level = "full_sib"} or a
#' half-sibling level, where the among-group variance is generated purely by
#' family structure. Since mateR2 controls that structure directly through N_P,
#' SR and mean mates, the resulting F-statistics can be read against a known
#' truth rather than estimated blind.
#'
#' \strong{Allele encoding, and a trap worth knowing about.} hierfstat infers the
#' number of digits per allele from the largest genotype value in the whole
#' table, not per locus (see \code{hierfstat:::getal.b}). A one-digit locus
#' therefore decodes correctly on its own but silently mis-decodes as soon as any
#' allele anywhere reaches 10. This function sidesteps that by choosing one width
#' for the entire table: two digits when the largest allele is at most 99, three
#' when it is at most 999. The resulting maximum always lands in the range that
#' forces hierfstat to the matching modulo. Alleles are written in ascending
#' order within a genotype, which keeps hierfstat's three-digit heuristic on the
#' correct branch.
#'
#' \strong{Missing genotypes.} COLONY writes a missing allele as \code{0}.
#' Encoding that would produce a genotype indistinguishable from a real one, so
#' any locus with a zero allele becomes \code{NA}, which is what hierfstat reads
#' as missing.
#'
#' \strong{Dominant markers} carry no heterozygote information and are refused
#' rather than silently treated as codominant.
#'
#' @param dat An object from \code{\link{read_colony_dat}}.
#' @param level Grouping for the first column. One of \code{"full_sib"}
#'   (dam and sire pair), \code{"paternal"}, \code{"maternal"}, or
#'   \code{"cohort"}.
#' @param min_size Drop groups with fewer than this many offspring. Within-group
#'   estimators are unstable at very small n, and a family of one contributes no
#'   within-family variance at all. Default 1 keeps everything.
#' @param labels Logical. Attach the group labels as a \code{"pop_labels"}
#'   attribute mapping the integer codes back to family identifiers.
#'
#' @return A data frame whose first column \code{pop} is an integer group code
#'   and whose remaining columns are one per locus, ready for
#'   \code{hierfstat::basic.stats()}, \code{wc()} or \code{genet.dist()}.
#'   Attributes \code{pop_labels} and \code{ncode} record the group mapping and
#'   the digits per allele.
#'
#' @seealso \code{\link{colony_hierfstat_levels}} for the nested design that
#'   \code{hierfstat::varcomp.glob()} takes.
#'
#' @examples
#' \dontrun{
#' dat <- read_colony_dat("COLONY2_1.DAT", id_key = "Input3_id_key.csv")
#' hf  <- colony_to_hierfstat(dat, level = "full_sib", min_size = 2)
#' hierfstat::basic.stats(hf)$overall
#' }
#'
#' @export
colony_to_hierfstat <- function(dat,
                                level = c("full_sib", "paternal", "maternal", "cohort"),
                                min_size = 1,
                                labels = TRUE) {

  level <- match.arg(level)
  if (!inherits(dat, "mateR2_colony_dat")) {
    stop("\n[mateR2 HIERFSTAT EXPORT ERROR]\n'dat' must come from read_colony_dat().",
         call. = FALSE)
  }
  if (any(dat$markers$type == 1)) {
    stop(sprintf(
      "\n[mateR2 HIERFSTAT EXPORT ERROR]\n%d of %d loci are dominant. hierfstat's estimators need codominant genotypes,\nbecause a dominant marker cannot distinguish a heterozygote from a dominant homozygote.\nSimulate with marker_type = 0 to use these statistics.",
      sum(dat$markers$type == 1), nrow(dat$markers)), call. = FALSE)
  }

  grp <- .hf_group(dat, level)
  keep <- rep(TRUE, length(grp))
  if (min_size > 1) {
    tab  <- table(grp)
    keep <- grp %in% names(tab)[tab >= min_size]
    if (!any(keep)) {
      stop(sprintf(
        "\n[mateR2 HIERFSTAT EXPORT ERROR]\nNo %s group reaches min_size = %d; the largest holds %d offspring.\nLower min_size, or sample more offspring per family.",
        level, min_size, max(tab)), call. = FALSE)
    }
  }

  gt  <- as.matrix(dat$offspring[keep, -1, drop = FALSE])
  grp <- grp[keep]
  enc <- .hf_encode(gt, nrow(dat$markers))

  pop_f <- factor(grp, levels = unique(grp))
  out <- data.frame(pop = as.integer(pop_f))
  out <- cbind(out, as.data.frame(enc$geno))
  names(out)[-1] <- dat$markers$marker
  rownames(out) <- NULL

  if (labels) attr(out, "pop_labels") <- levels(pop_f)
  attr(out, "ncode") <- enc$ncode
  attr(out, "level") <- level
  out
}


#' @title Nested Levels for a Hierarchical Variance Decomposition
#'
#' @description Builds the levels data frame that
#'   \code{hierfstat::varcomp.glob()} takes, nesting offspring within full-sibling
#'   family within one parental half-sibling family.
#'
#' @details
#' \strong{Only one sex can be the outer level.} A full-sibling family has
#' exactly one sire and exactly one dam, so nesting full-sib families inside
#' paternal families is well defined, and so is nesting them inside maternal
#' families. Nesting inside both at once is not: as soon as the mating network
#' has any multiple mating, paternal and maternal families cross rather than
#' nest, and a sire's families are spread across several dams. That is the whole
#' point of the bipartite structure mateR2 simulates, and it is why this function
#' takes \code{outer} rather than returning both. Running the same cohort twice,
#' once with each sex outer, is the informative comparison, and the difference
#' between the two is a direct read on the sex asymmetry in the mating system.
#'
#' @param dat An object from \code{\link{read_colony_dat}}.
#' @param outer Which parental sex forms the outer level, \code{"paternal"}
#'   (default) or \code{"maternal"}.
#' @param min_size Drop full-sibling families smaller than this, matched to the
#'   \code{min_size} passed to \code{\link{colony_to_hierfstat}} so the rows line
#'   up.
#'
#' @return A data frame of two factor columns, outermost first: the parental
#'   family and the full-sibling family, with one row per retained offspring.
#'
#' @examples
#' \dontrun{
#' hf  <- colony_to_hierfstat(dat, level = "full_sib", min_size = 2)
#' lev <- colony_hierfstat_levels(dat, outer = "paternal", min_size = 2)
#' hierfstat::varcomp.glob(levels = lev, loci = hf[, -1])$F
#' }
#'
#' @export
colony_hierfstat_levels <- function(dat, outer = c("paternal", "maternal"),
                                    min_size = 1) {
  outer <- match.arg(outer)
  if (!inherits(dat, "mateR2_colony_dat")) {
    stop("\n[mateR2 HIERFSTAT EXPORT ERROR]\n'dat' must come from read_colony_dat().",
         call. = FALSE)
  }
  fs <- .hf_group(dat, "full_sib")
  keep <- rep(TRUE, length(fs))
  if (min_size > 1) {
    tab  <- table(fs)
    keep <- fs %in% names(tab)[tab >= min_size]
  }
  par <- .hf_group(dat, outer)
  out <- data.frame(parent   = factor(par[keep]),
                    full_sib = factor(fs[keep]),
                    stringsAsFactors = FALSE)

  # A parent that mated once contributes a single full-sibling family, so the
  # intermediate level has no variance to estimate. varcomp.glob() reports NaN
  # rather than complaining, so catch it here where the cause is visible.
  per_parent <- tapply(as.character(out$full_sib), out$parent,
                       function(v) length(unique(v)))
  other <- if (outer == "paternal") "maternal" else "paternal"
  if (all(per_parent == 1L)) {
    stop(sprintf(
      "\n[mateR2 HIERFSTAT EXPORT ERROR]\nEvery %s family holds exactly one full-sibling family, so full-sib families are not nested within them and the intermediate level has no variance to estimate.\nThis happens when the sex forming the outer level mated only once. Try outer = \"%s\", or simulate a network in which %s parents have more than one mate.",
      outer, other, outer), call. = FALSE)
  }
  if (mean(per_parent == 1L) > 0.5) {
    warning(sprintf(
      "%d of %d %s families hold a single full-sibling family, so the intermediate level is weakly identified. Consider outer = \"%s\".",
      sum(per_parent == 1L), length(per_parent), outer, other), call. = FALSE)
  }
  out
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.hf_group <- function(dat, level) {
  t <- dat$truth
  dad <- if (!is.null(t$dad) && !all(is.na(t$dad))) t$dad else paste0("M", t$dad_index)
  mum <- if (!is.null(t$mum) && !all(is.na(t$mum))) t$mum else paste0("F", t$mum_index)
  switch(level,
         full_sib = paste(dad, mum, sep = "x"),
         paternal = dad,
         maternal = mum,
         cohort   = rep("cohort", nrow(t)))
}

# Encode a wide allele matrix (two columns per locus) into hierfstat's single
# number per genotype. The digit width is fixed across the whole table on
# purpose; see the note in the roxygen about getal.b reading the global maximum.
.hf_encode <- function(gt, n_loci) {
  a1 <- gt[, seq(1, ncol(gt), by = 2), drop = FALSE]
  a2 <- gt[, seq(2, ncol(gt), by = 2), drop = FALSE]

  max_allele <- suppressWarnings(max(c(a1, a2), na.rm = TRUE))
  if (!is.finite(max_allele) || max_allele > 999) {
    stop(sprintf(
      "\n[mateR2 HIERFSTAT EXPORT ERROR]\nThe largest allele is %s. hierfstat encodes at most three digits per allele.",
      format(max_allele)), call. = FALSE)
  }
  ncode  <- if (max_allele <= 99) 2L else 3L
  modulo <- if (ncode == 2L) 100L else 1000L

  # A zero in either position is COLONY's missing code; the genotype is missing.
  miss <- (a1 == 0) | (a2 == 0) | is.na(a1) | is.na(a2)

  lo <- pmin(a1, a2)
  hi <- pmax(a1, a2)
  geno <- lo * modulo + hi
  geno[miss] <- NA_integer_
  storage.mode(geno) <- "integer"
  colnames(geno) <- NULL

  list(geno = geno, ncode = ncode)
}
