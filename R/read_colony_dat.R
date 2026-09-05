#' @title Read a COLONY2 Data File
#'
#' @description Parses the \code{Colony2.dat} file written by COLONY's
#'   simulation module (\code{Simu2.exe}) back into R: run parameters, marker
#'   definitions, offspring and candidate-parent genotypes, and the true pedigree.
#'   This is the return leg of \code{\link{write_colony_sim}} and closes the
#'   benchmarking loop, since the recovered truth is what a reconstruction is
#'   scored against.
#'
#' @details
#' \strong{Where the truth comes from.} The simulation module encodes the true
#' pedigree in the offspring identifiers it generates: \code{M6F4C3} is the third
#' child of the male in row 6 and the female in column 4 of the mating matrix
#' supplied in \code{Input3.Par}. Those indices are parsed into
#' \code{dad_index} and \code{mum_index}. Pass the key written by
#' \code{write_colony_sim(write_key = TRUE)} as \code{id_key} and the mateR2
#' identifiers are joined on, so the truth table reads \code{M001}/\code{F001}
#' rather than positions. The file's own \code{True Configuration} block is
#' returned separately as \code{config} and is a useful cross-check: it gives
#' paternal and maternal family assignments per offspring independently of the
#' identifier convention.
#'
#' \strong{Candidate parents.} COLONY distinguishes the sampled true parents from
#' the unrelated padding by case: \code{M1}, \code{M2}, ... are true parents
#' carried over from the mating matrix, \code{m7}, \code{m8}, ... are the
#' unrelated candidates. That is recorded in the \code{true_parent} column rather
#' than left for the caller to infer from capitalisation.
#'
#' \strong{Genotype layout.} Genotypes are returned wide, one row per individual,
#' with the identifier in the first column and two columns per locus named
#' \code{Mk1-1}, \code{Mk1-2} and so on, which is the layout COLONY itself uses
#' and the layout a \code{.dat} writer expects. Missing alleles are COLONY's
#' \code{0}, left as-is rather than converted to \code{NA}, because a
#' \code{.dat} writer needs the zeros back.
#'
#' @param path Path to the \code{.DAT} file, for example \code{COLONY2_1.DAT}.
#' @param id_key Optional. The identifier key written by
#'   \code{\link{write_colony_sim}}, either as a path to the CSV or the data
#'   frame itself. When supplied, mateR2 identifiers are joined onto the truth
#'   table and the candidate-parent tables.
#'
#' @return A list of class \code{mateR2_colony_dat} with elements:
#'   \describe{
#'     \item{params}{Named list of the run parameters from the file header.}
#'     \item{markers}{Data frame with one row per locus: \code{marker},
#'       \code{type}, \code{dropout}, \code{error}.}
#'     \item{offspring}{Offspring genotypes, wide.}
#'     \item{candidate_males, candidate_females}{Candidate-parent genotypes,
#'       wide, with a \code{true_parent} flag. Zero rows under sibship-only.}
#'     \item{truth}{One row per offspring: \code{off}, \code{dad_index},
#'       \code{mum_index}, \code{sib_index}, plus \code{dad}/\code{mum} mateR2
#'       identifiers when \code{id_key} is supplied.}
#'     \item{config}{The file's own True Configuration block.}
#'     \item{path}{The file read.}
#'   }
#'
#' @seealso \code{\link{write_colony_sim}} for the outward leg,
#'   \code{\link{colony_truth_matrix}} to rebuild the mating matrix from the
#'   recovered truth.
#'
#' @examples
#' \dontrun{
#' dat <- read_colony_dat("COLONY2_1.DAT", id_key = "Input3_id_key.csv")
#' head(dat$offspring[, 1:5])
#' head(dat$truth)
#' }
#'
#' @importFrom utils read.csv
#' @export
read_colony_dat <- function(path, id_key = NULL) {

  if (!file.exists(path)) {
    stop(sprintf("\n[mateR2 COLONY READ ERROR]\n'%s' does not exist.", path),
         call. = FALSE)
  }
  L <- sub("\r$", "", readLines(path, warn = FALSE))

  # --------------------------------------------------------------------------
  # 1. Header. Every parameter line carries a "!" note; the marker block is the
  #    first line without one. Scanning rather than counting keeps this working
  #    if COLONY's header ever gains or loses a field.
  # --------------------------------------------------------------------------
  is_param <- grepl("!", L, fixed = TRUE)
  last_param <- which(!is_param & nzchar(trimws(L)))[1] - 1L
  if (is.na(last_param) || last_param < 5L) {
    stop(sprintf("\n[mateR2 COLONY READ ERROR]\n'%s' does not look like a COLONY2 data file; no parameter header found.", path),
         call. = FALSE)
  }
  hdr <- L[seq_len(last_param)]
  params <- .colony_parse_params(hdr)

  # --------------------------------------------------------------------------
  # 2. Marker block: names, types, dropout rates, other-error rates.
  # --------------------------------------------------------------------------
  pos <- last_param + 1L
  mk <- lapply(0:3, function(k) .colony_split(L[pos + k]))
  n_loci <- length(mk[[1]])
  if (!all(lengths(mk) == n_loci)) {
    stop(sprintf("\n[mateR2 COLONY READ ERROR]\nThe marker block is ragged: %s values across the four lines.",
                 paste(lengths(mk), collapse = "/")), call. = FALSE)
  }
  markers <- data.frame(
    marker  = mk[[1]],
    type    = as.integer(mk[[2]]),
    dropout = as.numeric(mk[[3]]),
    error   = as.numeric(mk[[4]]),
    stringsAsFactors = FALSE
  )
  gt_names <- c(rbind(paste0(markers$marker, "-1"), paste0(markers$marker, "-2")))
  pos <- pos + 4L

  # --------------------------------------------------------------------------
  # 3. Offspring genotypes, running to the first blank line.
  # --------------------------------------------------------------------------
  blk <- .colony_take_block(L, pos)
  offspring <- .colony_genotype_frame(blk$lines, gt_names, "Offspring", path)
  pos <- blk$next_pos

  # --------------------------------------------------------------------------
  # 4. Candidate parents. The counts line tells us how many of each to read.
  # --------------------------------------------------------------------------
  pos <- .colony_skip_blank(L, pos)
  pos <- pos + 1L                                   # inclusion-probability line
  pos <- .colony_skip_blank(L, pos)
  counts <- as.integer(.colony_split(sub("!.*$", "", L[pos])))
  if (length(counts) != 2L || anyNA(counts)) {
    stop(sprintf("\n[mateR2 COLONY READ ERROR]\nExpected a candidate male/female count line at line %d of '%s', found: %s",
                 pos, path, sQuote(trimws(L[pos]))), call. = FALSE)
  }
  n_cand_m <- counts[1]; n_cand_f <- counts[2]
  pos <- pos + 1L

  cand_m <- .colony_empty_genotypes(gt_names, "Male")
  cand_f <- .colony_empty_genotypes(gt_names, "Female")
  if (n_cand_m > 0L) {
    pos    <- .colony_skip_blank(L, pos)
    blk    <- .colony_take_block(L, pos, n_cand_m)
    cand_m <- .colony_genotype_frame(blk$lines, gt_names, "Male", path)
    pos    <- blk$next_pos
  }
  if (n_cand_f > 0L) {
    pos    <- .colony_skip_blank(L, pos)
    blk    <- .colony_take_block(L, pos, n_cand_f)
    cand_f <- .colony_genotype_frame(blk$lines, gt_names, "Female", path)
  }
  # COLONY writes true parents capitalised and unrelated padding lower-case.
  cand_m$true_parent <- grepl("^M", cand_m$Male)
  cand_f$true_parent <- grepl("^F", cand_f$Female)

  # --------------------------------------------------------------------------
  # 5. True Configuration block, if the file carries one.
  # --------------------------------------------------------------------------
  cfg_at <- grep("True Configuration", L, fixed = TRUE)
  config <- NULL
  if (length(cfg_at) == 1L) {
    rows <- list()
    for (l in L[(cfg_at + 1L):length(L)]) {
      v <- suppressWarnings(as.integer(.colony_split(l)))
      if (length(v) != 6L || anyNA(v)) break
      rows[[length(rows) + 1L]] <- v
    }
    if (length(rows)) {
      m <- do.call(rbind, rows)
      config <- data.frame(offspring_index  = m[, 1],
                           paternal_family  = m[, 2],
                           maternal_family  = m[, 3],
                           cluster          = m[, 4],
                           stringsAsFactors = FALSE)
    }
  }

  # --------------------------------------------------------------------------
  # 6. Recover the true pedigree from the identifiers COLONY generated.
  # --------------------------------------------------------------------------
  truth <- .colony_truth_from_ids(offspring[[1]], path)
  if (!is.null(id_key)) {
    joined <- .colony_join_key(truth, cand_m, cand_f, id_key)
    truth  <- joined$truth
    cand_m <- joined$males
    cand_f <- joined$females
  }

  structure(list(params            = params,
                 markers           = markers,
                 offspring         = offspring,
                 candidate_males   = cand_m,
                 candidate_females = cand_f,
                 truth             = truth,
                 config            = config,
                 path              = path),
            class = "mateR2_colony_dat")
}


#' @title Rebuild a Mating Matrix from a Recovered COLONY Pedigree
#'
#' @description Turns the truth table returned by \code{\link{read_colony_dat}}
#'   back into a mating matrix of full-sibling counts, so it can be compared
#'   cell for cell against the \strong{K} that was exported.
#'
#' @param truth The \code{truth} element of a \code{\link{read_colony_dat}}
#'   result, or the whole result.
#' @param n_dads,n_mums Optional dimensions. Default to the largest index seen,
#'   which understates the matrix when a parent produced no offspring; supply
#'   \code{dim(K)} to get a frame that lines up with the original.
#'
#' @return An integer matrix, rows males and columns females.
#'
#' @examples
#' \dontrun{
#' dat <- read_colony_dat("COLONY2_1.DAT")
#' identical(colony_truth_matrix(dat, nrow(K), ncol(K)), K)
#' }
#'
#' @export
colony_truth_matrix <- function(truth, n_dads = NULL, n_mums = NULL) {
  if (inherits(truth, "mateR2_colony_dat")) truth <- truth$truth
  if (!is.data.frame(truth) ||
      !all(c("dad_index", "mum_index") %in% names(truth))) {
    stop("\n[mateR2 COLONY READ ERROR]\n'truth' must be the truth table from read_colony_dat().",
         call. = FALSE)
  }
  nd <- if (is.null(n_dads)) max(truth$dad_index) else n_dads
  nm <- if (is.null(n_mums)) max(truth$mum_index) else n_mums
  out <- matrix(0L, nd, nm)
  tab <- table(truth$dad_index, truth$mum_index)
  out[cbind(as.integer(rownames(tab))[row(tab)],
            as.integer(colnames(tab))[col(tab)])] <- as.integer(tab)
  out
}


#' @export
print.mateR2_colony_dat <- function(x, ...) {
  cat("COLONY2 data file:", basename(x$path), "\n")
  cat("  offspring        :", nrow(x$offspring), "\n")
  cat("  loci             :", nrow(x$markers),
      sprintf("(%d codominant, %d dominant)",
              sum(x$markers$type == 0), sum(x$markers$type == 1)), "\n")
  cat("  candidate males  :", nrow(x$candidate_males),
      sprintf("(%d true parents)", sum(x$candidate_males$true_parent)), "\n")
  cat("  candidate females:", nrow(x$candidate_females),
      sprintf("(%d true parents)", sum(x$candidate_females$true_parent)), "\n")
  cat("  true pedigree    :", length(unique(x$truth$dad_index)), "dads x",
      length(unique(x$truth$mum_index)), "mums,",
      nrow(unique(x$truth[, c("dad_index", "mum_index")])), "full-sib families\n")
  invisible(x)
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.colony_split <- function(x) {
  v <- trimws(unlist(strsplit(trimws(x), "[,\t ]+")))
  v[nzchar(v)]
}

.colony_skip_blank <- function(L, pos) {
  while (pos <= length(L) && !nzchar(trimws(L[pos]))) pos <- pos + 1L
  pos
}

# Take consecutive non-blank lines from pos, optionally capped at n.
.colony_take_block <- function(L, pos, n = NULL) {
  out <- character(0)
  while (pos <= length(L) && nzchar(trimws(L[pos]))) {
    out <- c(out, L[pos]); pos <- pos + 1L
    if (!is.null(n) && length(out) == n) break
  }
  list(lines = out, next_pos = pos)
}

.colony_genotype_frame <- function(lines, gt_names, id_col, path) {
  if (!length(lines)) return(.colony_empty_genotypes(gt_names, id_col))
  parts <- strsplit(lines, "[,\t ]+")
  wid   <- lengths(parts)
  if (length(unique(wid)) != 1L) {
    bad <- which(wid != wid[1])[1]
    stop(sprintf("\n[mateR2 COLONY READ ERROR]\nGenotype rows in '%s' have inconsistent widths; row %d has %d fields against %d elsewhere.",
                 path, bad, wid[bad], wid[1]), call. = FALSE)
  }
  if (wid[1] - 1L != length(gt_names)) {
    stop(sprintf("\n[mateR2 COLONY READ ERROR]\nGenotype rows in '%s' carry %d allele columns but the marker block declares %d loci (%d columns).",
                 path, wid[1] - 1L, length(gt_names) / 2, length(gt_names)),
         call. = FALSE)
  }
  m  <- do.call(rbind, parts)
  df <- data.frame(id = trimws(m[, 1]), stringsAsFactors = FALSE)
  gt <- suppressWarnings(matrix(as.integer(m[, -1, drop = FALSE]), nrow = nrow(m)))
  colnames(gt) <- gt_names
  df <- cbind(df, as.data.frame(gt, stringsAsFactors = FALSE))
  names(df)[1] <- id_col
  df
}

.colony_empty_genotypes <- function(gt_names, id_col) {
  df <- data.frame(id = character(0), stringsAsFactors = FALSE)
  for (nm in gt_names) df[[nm]] <- integer(0)
  names(df)[1] <- id_col
  df$true_parent <- logical(0)
  df
}

.colony_parse_params <- function(hdr) {
  val  <- trimws(sub("!.*$", "", hdr))
  note <- trimws(sub("^[^!]*!", "", hdr))
  # Give the fields the names the guide uses, falling back to the note text.
  key <- c("dataset_name", "output_name", "n_offspring", "n_loci", "seed",
           "update_allele_freq", "species", "inbreeding", "ploidy", "gamy",
           "clone", "scale_sibship", "sibship_prior", "known_allele_freq",
           "n_runs", "run_length", "monitor_by", "monitor_interval",
           "gui_mode", "likelihood", "precision")
  if (length(val) != length(key)) key <- paste0("field_", seq_along(val))
  out <- as.list(val); names(out) <- key
  attr(out, "notes") <- note
  # Numeric fields become numbers; the two names stay strings.
  for (i in seq_along(out)) {
    v <- suppressWarnings(as.numeric(.colony_split(out[[i]])))
    if (!anyNA(v) && length(v)) out[[i]] <- if (length(v) == 1L) v else v
  }
  out
}

# Offspring identifiers generated by Simu2.exe carry the true pedigree:
# M<dad row>F<mum column>C<child number>.
.colony_truth_from_ids <- function(ids, path) {
  m <- regmatches(ids, regexec("^[Mm]([0-9]+)[Ff]([0-9]+)[Cc]([0-9]+)$", ids))
  ok <- lengths(m) == 4L
  if (!any(ok)) {
    warning(sprintf("Offspring identifiers in '%s' do not follow COLONY's M<i>F<j>C<k> simulation convention, so the true pedigree could not be recovered.",
                    path), call. = FALSE)
    return(data.frame(off = ids, dad_index = NA_integer_, mum_index = NA_integer_,
                      sib_index = NA_integer_, stringsAsFactors = FALSE))
  }
  if (!all(ok)) {
    warning(sprintf("%d of %d offspring identifiers in '%s' do not follow the M<i>F<j>C<k> convention; their pedigree is NA.",
                    sum(!ok), length(ids), path), call. = FALSE)
  }
  grab <- function(k) vapply(seq_along(m), function(i)
    if (ok[i]) as.integer(m[[i]][k + 1L]) else NA_integer_, integer(1))
  data.frame(off = ids, dad_index = grab(1), mum_index = grab(2),
             sib_index = grab(3), stringsAsFactors = FALSE)
}

.colony_join_key <- function(truth, males, females, id_key) {
  if (is.character(id_key)) {
    if (!file.exists(id_key)) {
      stop(sprintf("\n[mateR2 COLONY READ ERROR]\nIdentifier key '%s' does not exist.", id_key),
           call. = FALSE)
    }
    id_key <- utils::read.csv(id_key, stringsAsFactors = FALSE)
  }
  need <- c("sex", "index", "mateR2_id")
  if (!is.data.frame(id_key) || !all(need %in% names(id_key))) {
    stop("\n[mateR2 COLONY READ ERROR]\n'id_key' must be the CSV written by write_colony_sim(write_key = TRUE), with columns sex, index and mateR2_id.",
         call. = FALSE)
  }
  km <- id_key[id_key$sex == "male", ]
  kf <- id_key[id_key$sex == "female", ]
  truth$dad <- km$mateR2_id[match(truth$dad_index, km$index)]
  truth$mum <- kf$mateR2_id[match(truth$mum_index, kf$index)]

  idx <- function(v) suppressWarnings(as.integer(sub("^[A-Za-z]+", "", v)))
  if (nrow(males)) {
    males$mateR2_id <- ifelse(males$true_parent,
                              km$mateR2_id[match(idx(males[[1]]), km$index)], NA)
  } else males$mateR2_id <- character(0)
  if (nrow(females)) {
    females$mateR2_id <- ifelse(females$true_parent,
                                kf$mateR2_id[match(idx(females[[1]]), kf$index)], NA)
  } else females$mateR2_id <- character(0)

  list(truth = truth, males = males, females = females)
}


#' @title Reshape a COLONY2 Data File into Inputs for a .dat Writer
#'
#' @description Turns the object returned by \code{\link{read_colony_dat}} into
#'   the four pieces a COLONY \code{.dat} writer takes: an offspring genotype
#'   frame, candidate father and mother frames, and a marker block. This is what
#'   lets a simulated cohort be handed straight to an inference run, optionally
#'   sub-sampled, with or without a candidate parent list.
#'
#' @details
#' The marker block is returned in the layout an inference \code{.dat} expects:
#' one column per locus named for the marker, and three rows giving marker type,
#' allelic dropout rate and other-error rate. The allele-count row that appears
#' in the simulation input is not part of an inference file and is dropped.
#'
#' \code{parents = "none"} produces empty parent frames, which is the sibship-only
#' case. \code{"true"} keeps only the real parents and discards the unrelated
#' padding, which is useful for measuring assignment accuracy against a clean
#' candidate list. \code{"all"} keeps the padding, which is the realistic case
#' and the one that exercises false-assignment rates.
#'
#' @param dat An object from \code{\link{read_colony_dat}}.
#' @param n_offspring Optional. Draw this many offspring at random without
#'   replacement, to emulate sampling effort. Defaults to all of them.
#' @param parents One of \code{"all"}, \code{"true"} or \code{"none"}.
#' @param use_mateR2_ids Logical. Label parents with their mateR2 identifiers
#'   where known, which requires that \code{read_colony_dat()} was given an
#'   \code{id_key}. Offspring keep COLONY's identifiers either way, since those
#'   carry the pedigree.
#'
#' @return A list with \code{kids}, \code{dads}, \code{moms}, \code{markers} and
#'   \code{truth}, the last subset to the offspring retained.
#'
#' @examples
#' \dontrun{
#' dat <- read_colony_dat("COLONY2_1.DAT", id_key = "Input3_id_key.csv")
#' x <- colony_dat_inputs(dat, n_offspring = 200, parents = "none")
#' str(x$kids[, 1:5])
#' }
#'
#' @export
colony_dat_inputs <- function(dat, n_offspring = NULL, parents = c("all", "true", "none"),
                              use_mateR2_ids = TRUE) {
  parents <- match.arg(parents)
  if (!inherits(dat, "mateR2_colony_dat")) {
    stop("\n[mateR2 COLONY READ ERROR]\n'dat' must come from read_colony_dat().",
         call. = FALSE)
  }

  kids  <- dat$offspring
  truth <- dat$truth
  if (!is.null(n_offspring)) {
    if (n_offspring > nrow(kids)) {
      stop(sprintf("\n[mateR2 COLONY READ ERROR]\nAsked for %d offspring but the file holds %d.",
                   n_offspring, nrow(kids)), call. = FALSE)
    }
    keep  <- sort(sample.int(nrow(kids), n_offspring))
    kids  <- kids[keep, , drop = FALSE]
    truth <- truth[keep, , drop = FALSE]
    rownames(kids) <- rownames(truth) <- NULL
  }

  trim <- function(df) {
    if (parents == "none") return(df[0, setdiff(names(df), c("true_parent", "mateR2_id")), drop = FALSE])
    if (parents == "true") df <- df[df$true_parent, , drop = FALSE]
    if (use_mateR2_ids && "mateR2_id" %in% names(df)) {
      df[[1]] <- ifelse(is.na(df$mateR2_id), df[[1]], df$mateR2_id)
    }
    df <- df[, setdiff(names(df), c("true_parent", "mateR2_id")), drop = FALSE]
    rownames(df) <- NULL
    df
  }

  # Inference .dat marker block: names as column headers, then type, dropout,
  # other-error. The simulation input's allele-count row has no place here.
  markers <- as.data.frame(rbind(dat$markers$type,
                                 dat$markers$dropout,
                                 dat$markers$error))
  names(markers) <- dat$markers$marker

  list(kids    = kids,
       dads    = trim(dat$candidate_males),
       moms    = trim(dat$candidate_females),
       markers = markers,
       truth   = truth)
}
