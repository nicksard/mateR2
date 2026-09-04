#' @title Write a COLONY Simulation-Module Input File (Input3.Par)
#'
#' @description Writes a weighted breeding matrix (\strong{K}) produced by
#'   \code{\link{brd.mat.fitness}} to an \code{Input3.Par} file for COLONY's
#'   simulation module (\code{Simu2.exe}; Wang 2013). COLONY's mating matrix is
#'   exactly the object mateR2 already produces: rows are male parents, columns
#'   are female parents, and cell \eqn{q_{ij}} is the number of full-sibling
#'   offspring produced by male \emph{i} with female \emph{j}. No restructuring
#'   is required.
#'
#' @details
#' \strong{Sibship-only versus parentage.} \code{mode = "sibship"} emits a file
#' with no candidate parents: every entry of the sampled-dads and sampled-mums
#' rows is 0 and the unrelated-candidate counts are forced to 0, so COLONY is
#' asked to reconstruct sibship from the offspring genotypes alone.
#' \code{mode = "parentage"} marks all true parents as sampled by default and
#' lets \code{n_candidate_males} / \code{n_candidate_females} pad the candidate
#' lists with unrelated individuals. Either default can be overridden per
#' individual through \code{sampled_dads} and \code{sampled_mums}.
#'
#' \strong{Identity.} COLONY generates its own individual identifiers and gives
#' the user no way to supply parent names, so the mateR2 identifiers carried on
#' the dimnames of \code{K} cannot travel into the file. Positional order is the
#' only link. With \code{write_key = TRUE} the function writes a companion CSV
#' mapping matrix row and column position back to the mateR2 identifier, which
#' is what makes scoring a reconstruction against ground truth possible.
#'
#' \strong{Writing.} The whole file is assembled in memory and committed with a
#' single \code{\link{writeLines}} call on one connection. Building the file in
#' many appending writes to the same path is what produces the stray blank lines
#' and truncated records that this format is prone to.
#'
#' \strong{Two contradictions in the COLONY user guide}, resolved here against
#' the guide's own worked example, which is self-consistent and matches observed
#' program behaviour:
#' \itemize{
#'   \item Section 6.1 item (12) reads as though 0 marks a sampled parent. The
#'     example annotates the same line "included (1) in or excluded (0) from
#'     candidate list", and section 5.3 item (12) confirms that a checked (that
#'     is, sampled) parent is included. \strong{1 means sampled.}
#'   \item Section 6.1 item (37) says to set the Windows/DOS indicator to 1 for a
#'     DOS run, but the example annotates it "0/1 for Windows DOS/GUI run" and
#'     uses 0. \strong{0 means DOS}, which is \code{gui_mode = FALSE} here.
#' }
#'
#' \strong{Verified against Simu2.exe.} The emitted file has been run through
#' COLONY's simulation module. On a 6 x 4 matrix in sibship mode the program
#' reports 40 offspring, 6 dads and 4 mums, writes \code{Colony2.dat} with
#' \code{0 0} candidate parents, and its own true-configuration block reproduces
#' the input matrix cell for cell. In parentage mode the candidate lists come
#' back as the true parents plus the requested unrelated padding.
#'
#' @param K Weighted breeding matrix, or anything coercible by
#'   \code{as.matrix()}. Rows are males, columns are females, cells are
#'   full-sibling counts. Whole numbers, non-negative. All-zero rows and columns
#'   are permitted and denote a parent contributing no offspring to the sample.
#' @param path File to write. COLONY's simulation module requires the name
#'   \code{Input3.Par} in the project folder (default).
#' @param mode Either \code{"sibship"} (default) or \code{"parentage"}. See
#'   Details.
#' @param project Output-file name stem used by COLONY, fewer than 40
#'   alphanumeric characters.
#' @param n_matrices Number of replicate mating matrices. COLONY simulates the
#'   structure this many times over distinct individuals, so total offspring is
#'   \code{n_matrices * sum(K)}. Leave at 1 to simulate the realized network.
#' @param n_replicates Number of replicate datasets to simulate and analyse.
#' @param n_loci Number of loci to simulate.
#' @param n_alleles Alleles per locus. Length 1 (recycled) or \code{n_loci}.
#' @param allele_dist One of \code{"uniform"}, \code{"equal"},
#'   \code{"triangular"}, \code{"empirical"}.
#' @param allele_freqs Required when \code{allele_dist = "empirical"}. A list of
#'   \code{n_loci} numeric vectors, or a matrix/data frame of \code{n_loci} rows.
#'   Vector \emph{i} must hold \code{n_alleles[i]} frequencies summing to 1.
#' @param dropout_rate,error_rate Per-locus allelic dropout and other-error
#'   rates. Length 1 (recycled) or \code{n_loci}.
#' @param marker_type Per-locus 0 for codominant, 1 for dominant. Length 1
#'   (recycled) or \code{n_loci}.
#' @param missing_rate Probability that a single-locus genotype is missing.
#' @param n_candidate_males,n_candidate_females Unrelated individuals padding the
#'   candidate lists. Must be 0 when \code{mode = "sibship"}.
#' @param prob_dad,prob_mom Assumed probability that a true father or mother is
#'   in the candidate list. Required by the format but only used when the
#'   candidate lists are non-empty.
#' @param sampled_dads,sampled_mums Optional 0/1 vectors of length
#'   \code{nrow(K)} and \code{ncol(K)} naming which true parents are genotyped
#'   and placed in the candidate lists. Defaults follow \code{mode}.
#' @param known_both,known_dad,known_mum Optional matrices matching
#'   \code{dim(K)} giving the number of offspring per cell with known paternity
#'   and maternity, known paternity only, and known maternity only. Default to
#'   all zeros. Constrained by
#'   \code{known_both + known_dad + known_mum <= K} cell-wise.
#' @param ploidy 1 for haplodiploid, 2 for diploid.
#' @param male_monogamous,female_monogamous Logical. \code{FALSE} (default) is
#'   polygamous, which is what a mateR2 network with mean mates above 1 requires.
#' @param selfing_rate Parental selfing rate, used for monoecious species only
#'   but required by the format.
#' @param sibship_prior 0/1/2/3 for no/weak/medium/strong sibship prior.
#' @param paternal_sib_size,maternal_sib_size Average sibship sizes accompanying
#'   the prior.
#' @param likelihood 0/1/2 for pairwise-likelihood score, full likelihood, or
#'   both.
#' @param precision 0/1/2/3 for low/medium/high/very high likelihood precision.
#' @param n_runs Replicate COLONY runs per dataset. The guide advises 1 for
#'   simulation.
#' @param run_length 0/1/2/3/4 for very short to very long.
#' @param map_length Genetic map length in Morgans; negative for infinite.
#' @param inbreeding_coef Inbreeding coefficient of parents.
#' @param inbreeding Logical. Allow inbreeding during COLONY inference.
#' @param clone Logical. Infer clones.
#' @param scale_sibship Logical. Scale full sibship.
#' @param known_allele_freqs Logical. Treat the simulating allele frequencies as
#'   known during inference.
#' @param update_allele_freqs Logical. Update allele frequencies from
#'   reconstructed pedigrees. Only meaningful when
#'   \code{known_allele_freqs = FALSE}.
#' @param gui_mode Logical. \code{FALSE} (default) writes the DOS indicator.
#' @param monitor_interval Iterations between progress reports.
#' @param seed Integer seed COLONY uses internally. When \code{NULL} (default)
#'   one is drawn from R's generator, so an upstream \code{set.seed()} fixes it.
#' @param wrap_at Optional integer. Split matrix rows across this many values per
#'   line. Matrix rows carry no trailing comment, so the format permits the
#'   split. Not normally needed: \code{Simu2.exe} was measured to read a
#'   6,002-character row from a 4 x 3000 matrix without complaint, and the widest
#'   row a mateR2 network produces at N_P = 1600 is roughly 1,100 characters.
#'   Kept as an escape hatch for wider matrices than that.
#' @param write_key Logical. Also write \code{<path>_id_key.csv} mapping matrix
#'   position to mateR2 identifier.
#' @param overwrite Logical. Overwrite \code{path} if it exists.
#'
#' @return Invisibly, a list with \code{path}, \code{key_path}, \code{lines}
#'   (the exact character vector written), and \code{n_offspring}.
#'
#' @references Wang, J. (2013) A simulation module in the computer program
#'   COLONY for sibship and parentage analysis. \emph{Molecular Ecology
#'   Resources} 13, 734-739.
#'
#' @seealso \code{\link{brd.mat.fitness}} for producing \code{K},
#'   \code{\link{mat.sub.sample}} for the sampled cohort.
#'
#' @examples
#' set.seed(1)
#' M <- matrix(0L, 6, 4,
#'             dimnames = list(sprintf("M%03d", 1:6), sprintf("F%03d", 1:4)))
#' M[cbind(1:6, c(1, 1, 2, 2, 3, 4))] <- 1L
#' K <- brd.mat.fitness(M, min.fert = 20, max.fert = 60)
#'
#' out <- write_colony_sim(K, path = tempfile(fileext = ".Par"),
#'                         n_loci = 12, n_alleles = 8)
#' head(out$lines, 8)
#'
#' @importFrom utils write.csv
#' @export
write_colony_sim <- function(K,
                             path = "Input3.Par",
                             mode = c("sibship", "parentage"),
                             project = "sim",
                             n_matrices = 1,
                             n_replicates = 1,
                             n_loci = 10,
                             n_alleles = 10,
                             allele_dist = c("uniform", "equal",
                                             "triangular", "empirical"),
                             allele_freqs = NULL,
                             dropout_rate = 0,
                             error_rate = 1e-04,
                             marker_type = 0,
                             missing_rate = 0,
                             n_candidate_males = 0,
                             n_candidate_females = 0,
                             prob_dad = 0.5,
                             prob_mom = 0.5,
                             sampled_dads = NULL,
                             sampled_mums = NULL,
                             known_both = NULL,
                             known_dad = NULL,
                             known_mum = NULL,
                             ploidy = 2,
                             male_monogamous = FALSE,
                             female_monogamous = FALSE,
                             selfing_rate = 0,
                             sibship_prior = 0,
                             paternal_sib_size = 1,
                             maternal_sib_size = 1,
                             likelihood = 1,
                             precision = 2,
                             n_runs = 1,
                             run_length = 2,
                             map_length = -1,
                             inbreeding_coef = 0,
                             inbreeding = FALSE,
                             clone = FALSE,
                             scale_sibship = TRUE,
                             known_allele_freqs = FALSE,
                             update_allele_freqs = FALSE,
                             gui_mode = FALSE,
                             monitor_interval = 1e+05,
                             seed = NULL,
                             wrap_at = NULL,
                             write_key = TRUE,
                             overwrite = FALSE) {

  mode        <- match.arg(mode)
  allele_dist <- match.arg(allele_dist)

  # --------------------------------------------------------------------------
  # 1. Validate the mating matrix
  # --------------------------------------------------------------------------
  K <- .colony_as_count_matrix(K, "K")
  n_dads <- nrow(K)
  n_mums <- ncol(K)

  if (n_dads < 1L || n_mums < 1L) {
    stop("\n[mateR2 COLONY EXPORT ERROR]\nThe mating matrix must have at least one row (male) and one column (female).",
         call. = FALSE)
  }
  if (sum(K) == 0) {
    stop("\n[mateR2 COLONY EXPORT ERROR]\nThe mating matrix contains no offspring; there is nothing for COLONY to simulate.",
         call. = FALSE)
  }

  # Known-parentage matrices: all zeros unless supplied.
  known_both <- .colony_known_matrix(known_both, dim(K), "known_both")
  known_dad  <- .colony_known_matrix(known_dad,  dim(K), "known_dad")
  known_mum  <- .colony_known_matrix(known_mum,  dim(K), "known_mum")

  # The guide constrains these cumulatively: r <= q, s + r <= q, t + s + r <= q.
  if (any(known_both + known_dad + known_mum > K)) {
    bad <- which(known_both + known_dad + known_mum > K, arr.ind = TRUE)[1, ]
    stop(sprintf(
      "\n[mateR2 COLONY EXPORT ERROR]\nThe known-parentage matrices claim more offspring than the mating matrix holds.\nAt row %d, column %d: known_both + known_dad + known_mum = %g but K = %g.\nCOLONY requires the known counts to sum to no more than the cell total.",
      bad[["row"]], bad[["col"]],
      known_both[bad[["row"]], bad[["col"]]] +
        known_dad[bad[["row"]], bad[["col"]]] +
        known_mum[bad[["row"]], bad[["col"]]],
      K[bad[["row"]], bad[["col"]]]
    ), call. = FALSE)
  }

  # --------------------------------------------------------------------------
  # 2. Validate the project name and per-locus vectors
  # --------------------------------------------------------------------------
  if (!is.character(project) || length(project) != 1L || is.na(project)) {
    stop("\n[mateR2 COLONY EXPORT ERROR]\n'project' must be a single character string.",
         call. = FALSE)
  }
  if (nchar(project) >= 40L || grepl("[^A-Za-z0-9]", project)) {
    stop(sprintf(
      "\n[mateR2 COLONY EXPORT ERROR]\nProject name '%s' is not acceptable to COLONY.\nIt must be fewer than 40 characters and contain only letters and numbers.",
      project
    ), call. = FALSE)
  }

  n_loci <- .colony_scalar_count(n_loci, "n_loci")
  n_alleles    <- .colony_recycle(n_alleles,    n_loci, "n_alleles")
  dropout_rate <- .colony_recycle(dropout_rate, n_loci, "dropout_rate")
  error_rate   <- .colony_recycle(error_rate,   n_loci, "error_rate")
  marker_type  <- .colony_recycle(marker_type,  n_loci, "marker_type")

  if (any(n_alleles < 2)) {
    stop("\n[mateR2 COLONY EXPORT ERROR]\nEvery locus needs at least 2 alleles.",
         call. = FALSE)
  }
  if (!all(marker_type %in% c(0, 1))) {
    stop("\n[mateR2 COLONY EXPORT ERROR]\n'marker_type' must be 0 (codominant) or 1 (dominant) at every locus.",
         call. = FALSE)
  }
  if (any(marker_type == 1 & n_alleles != 2)) {
    stop("\n[mateR2 COLONY EXPORT ERROR]\nDominant markers are fixed at 2 alleles; 'n_alleles' disagrees at one or more loci.",
         call. = FALSE)
  }
  for (nm in c("dropout_rate", "error_rate")) {
    v <- get(nm)
    if (any(v < 0 | v >= 1)) {
      stop(sprintf("\n[mateR2 COLONY EXPORT ERROR]\n'%s' must be at least 0 and less than 1 at every locus.", nm),
           call. = FALSE)
    }
  }
  if (missing_rate < 0 || missing_rate >= 1) {
    stop("\n[mateR2 COLONY EXPORT ERROR]\n'missing_rate' must be at least 0 and less than 1.",
         call. = FALSE)
  }

  # --------------------------------------------------------------------------
  # 3. Resolve sampling design from 'mode'
  # --------------------------------------------------------------------------
  n_candidate_males   <- .colony_scalar_count(n_candidate_males,   "n_candidate_males",   allow_zero = TRUE)
  n_candidate_females <- .colony_scalar_count(n_candidate_females, "n_candidate_females", allow_zero = TRUE)

  if (mode == "sibship") {
    if (n_candidate_males > 0 || n_candidate_females > 0) {
      stop("\n[mateR2 COLONY EXPORT ERROR]\nmode = \"sibship\" reconstructs from offspring alone, so the candidate parent lists must be empty.\nSet n_candidate_males and n_candidate_females to 0, or use mode = \"parentage\".",
           call. = FALSE)
    }
    default_flag <- 0L
  } else {
    default_flag <- 1L
  }

  sampled_dads <- .colony_flags(sampled_dads, n_dads, default_flag, "sampled_dads")
  sampled_mums <- .colony_flags(sampled_mums, n_mums, default_flag, "sampled_mums")

  if (mode == "sibship" && (any(sampled_dads == 1L) || any(sampled_mums == 1L))) {
    stop("\n[mateR2 COLONY EXPORT ERROR]\nmode = \"sibship\" cannot mark any true parent as sampled.\nDrop the sampled_dads / sampled_mums arguments, or use mode = \"parentage\".",
         call. = FALSE)
  }

  # A parent named in a known-parentage matrix is in the sample whether or not
  # it is flagged, so sibship-only and known parentage are incompatible.
  if (mode == "sibship" && any(known_both + known_dad + known_mum > 0)) {
    stop("\n[mateR2 COLONY EXPORT ERROR]\nKnown parentage places those parents in the sample, which contradicts mode = \"sibship\".\nUse mode = \"parentage\" if some parentage is known.",
         call. = FALSE)
  }

  # --------------------------------------------------------------------------
  # 4. Allele frequencies
  # --------------------------------------------------------------------------
  dist_code <- match(allele_dist,
                     c("uniform", "equal", "triangular", "empirical")) - 1L
  freq_lines <- character(0)
  if (allele_dist == "empirical") {
    freq_lines <- .colony_freq_lines(allele_freqs, n_loci, n_alleles)
  } else if (!is.null(allele_freqs)) {
    warning("'allele_freqs' is ignored unless allele_dist = \"empirical\".",
            call. = FALSE)
  }

  # --------------------------------------------------------------------------
  # 5. Assemble the file in memory
  # --------------------------------------------------------------------------
  if (is.null(seed)) seed <- sample.int(9999L, 1L)

  n_matrices   <- .colony_scalar_count(n_matrices,   "n_matrices")
  n_replicates <- .colony_scalar_count(n_replicates, "n_replicates")
  n_runs       <- .colony_scalar_count(n_runs,       "n_runs")

  lines <- c(
    # (1) - (7)
    .cl(project,                                   "Project (output file) name"),
    .cl(n_replicates,                              "Number of replicates"),
    .cl(likelihood,                                "0/1/2=PLS/FL/FL-PLS combined"),
    .cl(precision,                                 "0/1/2/3=low/medium/high/very high precision"),
    .cl(paste(2L, .num(selfing_rate)),             "2/1=Dioecious/Monoecious, selfing rate for monoecious"),
    .cl(n_matrices,                                "Number of mating matrices"),
    .cl(paste(n_dads, n_mums),                     "# dads & mums in a mating structure"),
    "",
    # (8) the mating matrix itself
    .colony_matrix_lines(K, wrap_at),
    "",
    # (9) - (11) known parentage, all zeros unless supplied
    .colony_matrix_lines(known_both, wrap_at),
    "",
    .colony_matrix_lines(known_dad, wrap_at),
    "",
    .colony_matrix_lines(known_mum, wrap_at),
    "",
    # (12) - (13) which true parents are genotyped
    .colony_int_row(sampled_dads, wrap_at),
    .colony_int_row(sampled_mums, wrap_at),
    "",
    # (14) - (17)
    .cl(paste(n_candidate_males, n_candidate_females),
        "# unrelated candidate males & females"),
    .cl(paste(.num(prob_dad), .num(prob_mom)),
        "Assumed prbs of fathers & mothers included in candidates"),
    .cl(n_loci,                                    "Number of loci"),
    .cl(.num(missing_rate),                        "Prob. of missing genotypes"),
    "",
    # (18) - (22)
    .cl(.numv(dropout_rate),                       "Drop rate for each locus"),
    .cl(.numv(error_rate),                         "OtherErrorRate for each locus"),
    .cl(.numv(marker_type),                        "Marker types, 0/1=codominant/dominant"),
    .cl(.numv(n_alleles),                          "# alleles at each locus"),
    "",
    .cl(dist_code, "0/1/2/3=Uniform/Equal/Triangular/Empirical allele freq. distr."),
    # (23) present only under an empirical distribution
    freq_lines,
    # (24) - (38)
    .cl(ploidy,                                    "1/n=HaploDiploid/n-ploid species"),
    .cl(paste(as.integer(male_monogamous), as.integer(female_monogamous)),
        "1/0=Monogamy/Polygamy for males & females"),
    .cl(seed,                                      "Seed for random number generator"),
    .cl(paste(sibship_prior, .num(paternal_sib_size), .num(maternal_sib_size)),
        "0/1/2/3=No/Weak/Medium/Strong sibship prior"),
    .cl(as.integer(clone),                         "0/1=Clone inference =No/Yes"),
    .cl(as.integer(scale_sibship),                 "0/1=Scale full sibship=No/Yes"),
    .cl(as.integer(known_allele_freqs),            "1/0 (Y/N) for known allele frequency"),
    .cl(as.integer(update_allele_freqs),           "1/0 for updating allele freq. or not"),
    .cl(n_runs,                                    "#replicate runs"),
    .cl(run_length,                                "0/1/2/3/4=VeryShort/Short/Medium/Long/VeryLong run"),
    .cl(.num(map_length),                          "Map length in Morgans. <0 for infinite length"),
    .cl(.num(inbreeding_coef),                     "Inbreeding coefficient of parents"),
    .cl(as.integer(inbreeding),                    "0/1=N/Y for allowing inbreeding in Colony"),
    .cl(as.integer(gui_mode),                      "0/1 for Windows DOS/GUI run"),
    .cl(format(monitor_interval, scientific = FALSE),
        "#iterates/#second for Windows DOS/GUI run")
  )

  # --------------------------------------------------------------------------
  # 6. Commit to disk in one write
  # --------------------------------------------------------------------------
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop(sprintf(
      "\n[mateR2 COLONY EXPORT ERROR]\n'%s' already exists.\nMove it, delete it, or call write_colony_sim(overwrite = TRUE).",
      path
    ), call. = FALSE)
  }

  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con)

  key_path <- NULL
  if (isTRUE(write_key)) {
    key_path <- .colony_write_key(K, path, sampled_dads, sampled_mums)
  }

  n_offspring <- sum(K) * n_matrices
  message(sprintf(
    "Wrote %s: %d dads x %d mums, %s offspring across %d mating matri%s (%s).",
    path, n_dads, n_mums, format(n_offspring, big.mark = ","),
    n_matrices, if (n_matrices == 1L) "x" else "ces", mode
  ))

  invisible(list(path        = path,
                 key_path    = key_path,
                 lines       = lines,
                 n_offspring = n_offspring))
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Append COLONY's trailing "!" note to a value. Notes are legal on any line that
# is not part of a wrapped block.
.cl <- function(value, note) paste(value, paste0("!", note), sep = "\t")

# Format a number without scientific notation or trailing zeros. COLONY reads
# free-format reals, but 1e-04 is not a Fortran-readable literal.
.num <- function(x) {
  # formatC pads a vector to a common width; trim it so rows carry no leading
  # whitespace and stay compact at N_P = 1600 scale.
  trimws(formatC(x, format = "fg", digits = 10, drop0trailing = TRUE))
}

.numv <- function(x) paste(.num(x), collapse = " ")

.colony_as_count_matrix <- function(x, nm) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) {
    stop(sprintf("\n[mateR2 COLONY EXPORT ERROR]\n'%s' must be a matrix or data frame.", nm),
         call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (anyNA(x)) {
    stop(sprintf("\n[mateR2 COLONY EXPORT ERROR]\n'%s' contains missing values; COLONY needs a complete matrix of counts.", nm),
         call. = FALSE)
  }
  if (any(x < 0)) {
    stop(sprintf("\n[mateR2 COLONY EXPORT ERROR]\n'%s' contains negative values; cells are offspring counts.", nm),
         call. = FALSE)
  }
  if (any(abs(x - round(x)) > 1e-8)) {
    stop(sprintf("\n[mateR2 COLONY EXPORT ERROR]\n'%s' contains fractional values; cells must be whole numbers of full siblings.\nRound the fecundity allocation before exporting.", nm),
         call. = FALSE)
  }
  round(x)
}

.colony_known_matrix <- function(x, dims, nm) {
  if (is.null(x)) return(matrix(0, dims[1], dims[2]))
  x <- .colony_as_count_matrix(x, nm)
  if (!identical(dim(x), dims)) {
    stop(sprintf(
      "\n[mateR2 COLONY EXPORT ERROR]\n'%s' is %d x %d but the mating matrix is %d x %d; they must match.",
      nm, nrow(x), ncol(x), dims[1], dims[2]
    ), call. = FALSE)
  }
  x
}

.colony_scalar_count <- function(x, nm, allow_zero = FALSE) {
  floor_val <- if (allow_zero) 0L else 1L
  if (length(x) != 1L || is.na(x) || !is.numeric(x) ||
      abs(x - round(x)) > 1e-8 || x < floor_val) {
    stop(sprintf(
      "\n[mateR2 COLONY EXPORT ERROR]\n'%s' must be a single whole number of at least %d.",
      nm, floor_val
    ), call. = FALSE)
  }
  as.integer(round(x))
}

.colony_recycle <- function(x, n, nm) {
  if (!is.numeric(x) || anyNA(x)) {
    stop(sprintf("\n[mateR2 COLONY EXPORT ERROR]\n'%s' must be numeric and complete.", nm),
         call. = FALSE)
  }
  if (length(x) == 1L) return(rep(x, n))
  if (length(x) != n) {
    stop(sprintf(
      "\n[mateR2 COLONY EXPORT ERROR]\n'%s' has %d values but there are %d loci.\nSupply one value to apply to every locus, or one value per locus.",
      nm, length(x), n
    ), call. = FALSE)
  }
  x
}

.colony_flags <- function(x, n, default, nm) {
  if (is.null(x)) return(rep(as.integer(default), n))
  if (is.logical(x)) x <- as.integer(x)
  if (length(x) != n || anyNA(x) || !all(x %in% c(0L, 1L))) {
    stop(sprintf(
      "\n[mateR2 COLONY EXPORT ERROR]\n'%s' must be %d values of 0 or 1 (1 = the parent is genotyped and placed in the candidate list).",
      nm, n
    ), call. = FALSE)
  }
  as.integer(x)
}

.colony_freq_lines <- function(freqs, n_loci, n_alleles) {
  if (is.null(freqs)) {
    stop("\n[mateR2 COLONY EXPORT ERROR]\nallele_dist = \"empirical\" requires 'allele_freqs'.\nSupply a list of one numeric vector per locus, or choose another distribution.",
         call. = FALSE)
  }
  if (is.matrix(freqs) || is.data.frame(freqs)) {
    freqs <- split(as.matrix(freqs), row(as.matrix(freqs)))
    freqs <- lapply(freqs, function(v) v[!is.na(v)])
  }
  if (!is.list(freqs) || length(freqs) != n_loci) {
    stop(sprintf(
      "\n[mateR2 COLONY EXPORT ERROR]\n'allele_freqs' must hold exactly %d loci; it holds %d.",
      n_loci, length(freqs)
    ), call. = FALSE)
  }
  for (i in seq_len(n_loci)) {
    v <- as.numeric(freqs[[i]])
    if (length(v) != n_alleles[i]) {
      stop(sprintf(
        "\n[mateR2 COLONY EXPORT ERROR]\nLocus %d has %d allele frequencies but n_alleles says %d.\nCOLONY requires the two to agree at every locus.",
        i, length(v), n_alleles[i]
      ), call. = FALSE)
    }
    if (any(v < 0 | v > 1)) {
      stop(sprintf("\n[mateR2 COLONY EXPORT ERROR]\nLocus %d has an allele frequency outside [0, 1].", i),
           call. = FALSE)
    }
    if (sum(v) < 0.99 || sum(v) > 1.01) {
      stop(sprintf(
        "\n[mateR2 COLONY EXPORT ERROR]\nAllele frequencies at locus %d sum to %.4f.\nCOLONY requires each locus to sum to between 0.99 and 1.01.",
        i, sum(v)
      ), call. = FALSE)
    }
  }
  vapply(freqs, function(v) paste(.num(as.numeric(v)), collapse = " "),
         character(1), USE.NAMES = FALSE)
}

# Rows of an integer matrix as whitespace-delimited lines. The paste over
# columns is vectorised across rows, which matters at N_P = 1600 where the
# matrix runs to hundreds of thousands of cells.
.colony_matrix_lines <- function(m, wrap_at = NULL) {
  if (is.null(wrap_at) && all(m == 0)) {
    # Every row is identical; build the string once.
    return(rep(paste(rep("0", ncol(m)), collapse = " "), nrow(m)))
  }
  if (is.null(wrap_at)) {
    return(do.call(paste, c(as.data.frame(m), list(sep = " "))))
  }
  unlist(lapply(seq_len(nrow(m)),
                function(i) .colony_int_row(m[i, ], wrap_at)),
         use.names = FALSE)
}

# A single vector as one line, or several if wrapping is requested. Wrapped
# blocks carry no "!" note, which the format forbids on continued content.
.colony_int_row <- function(v, wrap_at = NULL) {
  txt <- .num(v)
  if (is.null(wrap_at) || length(txt) <= wrap_at) {
    return(paste(txt, collapse = " "))
  }
  groups <- split(txt, ceiling(seq_along(txt) / wrap_at))
  vapply(groups, paste, character(1), collapse = " ", USE.NAMES = FALSE)
}

.colony_write_key <- function(K, path, sampled_dads, sampled_mums) {
  dads <- rownames(K)
  mums <- colnames(K)
  if (is.null(dads) || is.null(mums)) {
    warning("The mating matrix has no dimnames, so the identifier key records positions only.\nMatrices from mp_table_to_matrix() carry M001/F001 names through the pipeline.",
            call. = FALSE)
    dads <- sprintf("row%d", seq_len(nrow(K)))
    mums <- sprintf("col%d", seq_len(ncol(K)))
  }
  key <- data.frame(
    sex        = c(rep("male", nrow(K)), rep("female", ncol(K))),
    index      = c(seq_len(nrow(K)), seq_len(ncol(K))),
    mateR2_id  = c(dads, mums),
    offspring  = c(rowSums(K), colSums(K)),
    genotyped  = c(sampled_dads, sampled_mums),
    stringsAsFactors = FALSE
  )
  # Strip a trailing extension without pulling in the tools package; the
  # character class stops the match escaping into a dotted directory name.
  key_path <- paste0(sub("\\.[^./\\\\]*$", "", path), "_id_key.csv")
  utils::write.csv(key, key_path, row.names = FALSE)
  key_path
}
