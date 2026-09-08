#' @title Per-Locus Statistics and Allele Frequencies from a COLONY Data File
#'
#' @description Summarises the genotypes parsed by \code{\link{read_colony_dat}}
#'   one locus at a time: how many alleles were actually observed, their
#'   frequencies, observed and expected heterozygosity, and the rate of missing
#'   genotypes. Because a simulated cohort was generated from parameters the user
#'   supplied, the dropout and error rates requested are carried alongside, so
#'   what came back can be checked against what was asked for.
#'
#' @details
#' \strong{Offspring and parents are separate samples.} Candidate parents are an
#' independent draw from the same allele pool, so \code{scope = "both"} returns
#' both and a \code{group} column distinguishes them. Comparing the two is a
#' check on the simulation rather than a description of it: allele frequencies
#' that agree between generations indicate the simulator sampled as intended,
#' while a systematic difference points at the transmission step. Under
#' sibship-only export there are no candidate parents, and \code{"parents"} then
#' returns zero rows with a warning rather than failing.
#'
#' \strong{Expected heterozygosity} is the unbiased estimator,
#' \eqn{H_e = (n / (n - 1)) (1 - \sum p_i^2)}, over the \emph{n} genotyped
#' individuals at that locus. With one or no genotyped individual it is
#' \code{NA} rather than zero.
#'
#' \strong{Missing genotypes} are COLONY's \code{0} allele. An individual counts
#' as missing at a locus when either allele is zero, and missing individuals are
#' excluded from the frequency and heterozygosity calculations rather than being
#' counted as homozygotes.
#'
#' \strong{Scale.} The computation is vectorised across loci rather than looping,
#' so a marker panel in the thousands is not a problem for this function. The
#' practical ceiling sits upstream: COLONY is a full-likelihood method and the
#' wide genotype frame \code{\link{read_colony_dat}} returns carries two columns
#' per locus, which is where a very large panel becomes awkward in R.
#'
#' @param dat An object from \code{\link{read_colony_dat}}.
#' @param scope Which individuals to summarise: \code{"offspring"} (default),
#'   \code{"parents"}, or \code{"both"}.
#' @param freqs Logical. Return the allele-frequency table. Set \code{FALSE} for
#'   the per-locus summary alone.
#'
#' @return A list with:
#'   \describe{
#'     \item{loci}{One row per locus per group: \code{group}, \code{marker},
#'       \code{type}, \code{n_typed}, \code{n_missing}, \code{missing_rate},
#'       \code{n_alleles}, \code{Ho}, \code{He}, and the \code{dropout} and
#'       \code{error} rates the simulation was given.}
#'     \item{freqs}{Long: \code{group}, \code{marker}, \code{allele},
#'       \code{count}, \code{freq}. Present only when \code{freqs = TRUE}.}
#'   }
#'
#' @seealso \code{\link{read_colony_dat}} for the parser,
#'   \code{\link{colony_to_hierfstat}} for F-statistics across pedigree groups.
#'
#' @examples
#' \dontrun{
#' dat <- read_colony_dat("COLONY2_1.DAT")
#' ls  <- colony_locus_stats(dat, scope = "both")
#' head(ls$loci)
#'
#' # Did a uniform allele-frequency simulation come back uniform?
#' with(subset(ls$freqs, group == "offspring"), tapply(freq, marker, range))
#' }
#'
#' @export
colony_locus_stats <- function(dat, scope = c("offspring", "parents", "both"),
                               freqs = TRUE) {

  scope <- match.arg(scope)
  if (!inherits(dat, "mateR2_colony_dat")) {
    stop("\n[mateR2 COLONY READ ERROR]\n'dat' must come from read_colony_dat().",
         call. = FALSE)
  }

  groups <- list()
  if (scope %in% c("offspring", "both")) {
    groups$offspring <- dat$offspring
  }
  if (scope %in% c("parents", "both")) {
    m <- dat$candidate_males
    f <- dat$candidate_females
    if (nrow(m) == 0L && nrow(f) == 0L) {
      warning("This file carries no candidate parents, which is what a sibship-only export produces. No parent rows returned.",
              call. = FALSE)
    }
    if (nrow(m)) groups$male   <- m
    if (nrow(f)) groups$female <- f
  }
  if (!length(groups)) {
    return(list(loci = .cls_empty_loci(), freqs = if (freqs) .cls_empty_freqs()))
  }

  out_loci <- list()
  out_frq  <- list()
  for (g in names(groups)) {
    res <- .cls_one(groups[[g]], dat$markers, g, freqs)
    out_loci[[g]] <- res$loci
    if (freqs) out_frq[[g]] <- res$freqs
  }

  loci <- do.call(rbind, out_loci)
  rownames(loci) <- NULL
  res <- list(loci = loci)
  if (freqs) {
    fr <- do.call(rbind, out_frq)
    rownames(fr) <- NULL
    res$freqs <- fr
  }
  res
}


# ---------------------------------------------------------------------------
# Internal
# ---------------------------------------------------------------------------

.cls_empty_loci <- function() {
  data.frame(group = character(0), marker = character(0), type = integer(0),
             n_typed = integer(0), n_missing = integer(0),
             missing_rate = numeric(0), n_alleles = integer(0),
             Ho = numeric(0), He = numeric(0),
             dropout = numeric(0), error = numeric(0),
             stringsAsFactors = FALSE)
}

.cls_empty_freqs <- function() {
  data.frame(group = character(0), marker = character(0), allele = integer(0),
             count = integer(0), freq = numeric(0), stringsAsFactors = FALSE)
}

# One group. Vectorised across loci: allele counting is a single tabulate() over
# a (locus, allele) index rather than a loop, which is what keeps a panel of
# thousands of markers cheap.
.cls_one <- function(df, markers, group, want_freqs) {

  gt <- as.matrix(df[, seq_len(2L * nrow(markers)) + 1L, drop = FALSE])
  storage.mode(gt) <- "integer"
  n_loci <- nrow(markers)
  n_ind  <- nrow(gt)

  a1 <- gt[, seq(1L, ncol(gt), by = 2L), drop = FALSE]
  a2 <- gt[, seq(2L, ncol(gt), by = 2L), drop = FALSE]

  # COLONY writes 0 for a missing allele; a genotype is missing if either is 0.
  miss <- (a1 == 0L) | (a2 == 0L) | is.na(a1) | is.na(a2)
  n_missing <- colSums(miss)
  n_typed   <- n_ind - n_missing

  het <- (a1 != a2) & !miss
  Ho  <- ifelse(n_typed > 0, colSums(het) / n_typed, NA_real_)

  # --- allele counts, all loci in one pass ---------------------------------
  max_allele <- max(c(a1, a2), na.rm = TRUE)
  if (!is.finite(max_allele) || max_allele < 1) max_allele <- 1L
  locus_of <- rep(seq_len(n_loci), each = n_ind)          # column-major order
  keep     <- !as.vector(miss)
  idx <- c((locus_of[keep] - 1L) * max_allele + as.vector(a1)[keep],
           (locus_of[keep] - 1L) * max_allele + as.vector(a2)[keep])
  counts <- tabulate(idx, nbins = n_loci * max_allele)
  cm <- matrix(counts, nrow = max_allele, ncol = n_loci)  # allele x locus

  tot <- colSums(cm)
  p   <- sweep(cm, 2, pmax(tot, 1L), "/")
  n_alleles <- colSums(cm > 0L)

  # Unbiased He over genotyped individuals; undefined below two.
  hom <- colSums(p^2)
  He  <- ifelse(n_typed > 1, (n_typed / (n_typed - 1)) * (1 - hom), NA_real_)

  loci <- data.frame(
    group        = group,
    marker       = markers$marker,
    type         = markers$type,
    n_typed      = as.integer(n_typed),
    n_missing    = as.integer(n_missing),
    missing_rate = n_missing / n_ind,
    n_alleles    = as.integer(n_alleles),
    Ho           = as.numeric(Ho),
    He           = as.numeric(He),
    dropout      = markers$dropout,
    error        = markers$error,
    stringsAsFactors = FALSE
  )

  if (!want_freqs) return(list(loci = loci))

  obs <- which(cm > 0L, arr.ind = TRUE)                   # rows: allele, locus
  freqs <- data.frame(
    group  = group,
    marker = markers$marker[obs[, "col"]],
    allele = as.integer(obs[, "row"]),
    count  = as.integer(cm[obs]),
    freq   = as.numeric(p[obs]),
    stringsAsFactors = FALSE
  )
  freqs <- freqs[order(match(freqs$marker, markers$marker), freqs$allele), ]
  rownames(freqs) <- NULL

  list(loci = loci, freqs = freqs)
}
