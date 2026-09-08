#' @title Read Genotypes from a COLONY2 Data File
#'
#' @description Parses the \code{Colony2.dat} file written by COLONY's
#'   simulation module (\code{Simu2.exe}) into R: the marker definitions and the
#'   genotypes of the offspring and of any candidate parents. This is the return
#'   leg of \code{\link{write_colony_sim}}, and its output feeds
#'   \code{\link{colony_locus_stats}}.
#'
#' @details
#' \strong{Genotype layout.} Genotypes are returned wide, one row per
#' individual, with the identifier in the first column and two columns per locus
#' named \code{Mk1-1}, \code{Mk1-2} and so on, which is the layout COLONY itself
#' uses. Missing alleles are left as COLONY's \code{0} rather than converted to
#' \code{NA}, so a writer can round-trip them;
#' \code{\link{colony_locus_stats}} treats a zero in either position as a
#' missing genotype.
#'
#' \strong{Candidate parents.} COLONY distinguishes sampled true parents from
#' unrelated padding by case, writing \code{M1} for the first and \code{m7} for
#' the second. Both are returned as given. A sibship-only file has no candidate
#' parents and both frames come back with zero rows.
#'
#' @param path Path to the \code{.DAT} file, for example \code{COLONY2_1.DAT}.
#'
#' @return A list of class \code{mateR2_colony_dat} with:
#'   \describe{
#'     \item{markers}{One row per locus: \code{marker}, \code{type},
#'       \code{dropout}, \code{error}.}
#'     \item{offspring}{Offspring genotypes, wide.}
#'     \item{candidate_males, candidate_females}{Candidate-parent genotypes,
#'       wide. Zero rows under sibship-only.}
#'     \item{path}{The file read.}
#'   }
#'
#' @seealso \code{\link{write_colony_sim}} for the outward leg,
#'   \code{\link{colony_locus_stats}} for per-locus summaries.
#'
#' @examples
#' \dontrun{
#' dat <- read_colony_dat("COLONY2_1.DAT")
#' head(dat$offspring[, 1:5])
#' colony_locus_stats(dat$offspring, dat$markers)
#' }
#'
#' @export
read_colony_dat <- function(path) {

  if (!file.exists(path)) {
    stop(sprintf("\n[mateR2 COLONY READ ERROR]\n'%s' does not exist.", path),
         call. = FALSE)
  }
  L <- sub("\r$", "", readLines(path, warn = FALSE))

  # Every parameter line in the header carries a "!" note; the marker block is
  # the first line without one. Scanning rather than counting keeps this working
  # if COLONY's header gains or loses a field.
  is_param   <- grepl("!", L, fixed = TRUE)
  last_param <- which(!is_param & nzchar(trimws(L)))[1] - 1L
  if (is.na(last_param) || last_param < 5L) {
    stop(sprintf("\n[mateR2 COLONY READ ERROR]\n'%s' does not look like a COLONY2 data file; no parameter header found.", path),
         call. = FALSE)
  }

  # --- Marker block: names, types, dropout rates, other-error rates ---------
  pos <- last_param + 1L
  mk  <- lapply(0:3, function(k) .colony_split(L[pos + k]))
  n_loci <- length(mk[[1]])
  if (!all(lengths(mk) == n_loci)) {
    stop(sprintf("\n[mateR2 COLONY READ ERROR]\nThe marker block is ragged: %s values across the four lines.",
                 paste(lengths(mk), collapse = "/")), call. = FALSE)
  }
  markers <- data.frame(marker  = mk[[1]],
                        type    = as.integer(mk[[2]]),
                        dropout = as.numeric(mk[[3]]),
                        error   = as.numeric(mk[[4]]),
                        stringsAsFactors = FALSE)
  gt_names <- c(rbind(paste0(markers$marker, "-1"), paste0(markers$marker, "-2")))
  pos <- pos + 4L

  # --- Offspring genotypes, to the first blank line -------------------------
  blk       <- .colony_take_block(L, pos)
  offspring <- .colony_genotype_frame(blk$lines, gt_names, "Offspring", path)
  pos       <- blk$next_pos

  # --- Candidate parents; the counts line says how many of each -------------
  pos <- .colony_skip_blank(L, pos)
  pos <- pos + 1L                                   # inclusion-probability line
  pos <- .colony_skip_blank(L, pos)
  counts <- as.integer(.colony_split(sub("!.*$", "", L[pos])))
  if (length(counts) != 2L || anyNA(counts)) {
    stop(sprintf("\n[mateR2 COLONY READ ERROR]\nExpected a candidate male/female count line at line %d of '%s', found: %s",
                 pos, path, sQuote(trimws(L[pos]))), call. = FALSE)
  }
  pos <- pos + 1L

  cand_m <- .colony_empty_genotypes(gt_names, "Male")
  cand_f <- .colony_empty_genotypes(gt_names, "Female")
  if (counts[1] > 0L) {
    pos    <- .colony_skip_blank(L, pos)
    blk    <- .colony_take_block(L, pos, counts[1])
    cand_m <- .colony_genotype_frame(blk$lines, gt_names, "Male", path)
    pos    <- blk$next_pos
  }
  if (counts[2] > 0L) {
    pos    <- .colony_skip_blank(L, pos)
    blk    <- .colony_take_block(L, pos, counts[2])
    cand_f <- .colony_genotype_frame(blk$lines, gt_names, "Female", path)
  }

  structure(list(markers           = markers,
                 offspring         = offspring,
                 candidate_males   = cand_m,
                 candidate_females = cand_f,
                 path              = path),
            class = "mateR2_colony_dat")
}


#' @export
print.mateR2_colony_dat <- function(x, ...) {
  cat("COLONY2 data file:", basename(x$path), "\n")
  cat("  offspring        :", nrow(x$offspring), "\n")
  cat("  loci             :", nrow(x$markers),
      sprintf("(%d codominant, %d dominant)",
              sum(x$markers$type == 0), sum(x$markers$type == 1)), "\n")
  cat("  candidate males  :", nrow(x$candidate_males), "\n")
  cat("  candidate females:", nrow(x$candidate_females), "\n")
  invisible(x)
}


# ---------------------------------------------------------------------------
# Internal
# ---------------------------------------------------------------------------

.colony_split <- function(x) {
  v <- trimws(unlist(strsplit(trimws(x), "[,\t ]+")))
  v[nzchar(v)]
}

.colony_skip_blank <- function(L, pos) {
  while (pos <= length(L) && !nzchar(trimws(L[pos]))) pos <- pos + 1L
  pos
}

# Consecutive non-blank lines from pos, optionally capped at n.
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
  gt <- suppressWarnings(matrix(as.integer(m[, -1, drop = FALSE]), nrow = nrow(m)))
  colnames(gt) <- gt_names
  df <- cbind(data.frame(id = trimws(m[, 1]), stringsAsFactors = FALSE),
              as.data.frame(gt, stringsAsFactors = FALSE))
  names(df)[1] <- id_col
  df
}

.colony_empty_genotypes <- function(gt_names, id_col) {
  df <- data.frame(id = character(0), stringsAsFactors = FALSE)
  for (nm in gt_names) df[[nm]] <- integer(0)
  names(df)[1] <- id_col
  df
}
