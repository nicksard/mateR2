#' @title Sub-Sample Offspring from a Breeding Matrix (Null Model)
#' @description This function serves as the Null Model for field sampling. It
#'   simulates the process of capturing a fixed number of offspring from a
#'   breeding matrix under the assumption of perfect panmixia (no
#'   spatial-temporal clustering). Offspring are drawn \strong{without
#'   replacement} from the realised cohort, so a mate pair can never contribute
#'   more sampled offspring than it actually produced. This represents the
#'   "best-case scenario" for sampling offspring for pedigree analysis and
#'   serves as a baseline against which empirical or clustered sampling methods
#'   can be tested.
#'
#' @details
#' The draw is an exact multivariate hypergeometric sample over the cells of the
#' breeding matrix: each mate pair's sampled count is drawn conditionally on the
#' offspring and sample budget still remaining. This is equivalent to shuffling
#' the whole cohort and taking the first `num_offspring` juveniles, but does not
#' require the cohort to be expanded in memory. Every offspring therefore has the
#' same inclusion probability, `num_offspring / sum(mat)`, and the expected number
#' drawn from a pair is proportional to that pair's fecundity.
#'
#' Mate pairs that produced offspring but had none sampled are retained in the
#' output with `off1 = 0`, which keeps undetected parents enumerable for
#' downstream yield-class statistics.
#'
#' @param mat The breeding matrix with offspring counts.
#' @param num_offspring The number of offspring to sub-sample. Cannot exceed the
#'   total cohort size `sum(mat)`.
#' @return A data frame in long format containing the sampled mate pairs,
#'   their original offspring count, and the new sampled count.
#' @importFrom stats rhyper
#' @export
mat.sub.sample <- function(mat, num_offspring) {
  # Validates inputs
  if (any(mat < 0) || any(mat %% 1 != 0)) {
    stop("Input matrix 'mat' must contain non-negative integers.")
  }
  if (num_offspring <= 0) {
    stop("'num_offspring' must be a positive integer.")
  }
  total_offspring <- sum(mat)
  if (num_offspring > total_offspring) {
    stop(sprintf(
      paste0("\n[mateR2 SAMPLING ERROR]\nRequested %.0f offspring but the cohort contains only %.0f.\n",
             "Sampling is without replacement, so the sample cannot exceed the cohort.\n",
             "Lower 'num_offspring', or raise fecundity when building the breeding matrix."),
      num_offspring, total_offspring
    ), call. = FALSE)
  }

  # Convert matrix to long format with mate pair IDs (using base R)
  ped1 <- as.data.frame(as.table(mat))
  colnames(ped1) <- c("dads", "moms", "off")
  ped1 <- ped1[ped1$off != 0, ]
  ped1$mp <- paste(ped1$dads, ped1$moms, sep = "_")

  # Per-pair sampling weight, retained for reference and downstream code
  ped1$probs <- ped1$off / total_offspring

  # --- Multivariate hypergeometric draw (sampling without replacement) ---
  # Walk the pairs once, drawing each pair's contribution conditional on the
  # offspring and sample budget still outstanding. rhyper(1, m, n, k) draws the
  # number of "white balls" among k draws from an urn of m white and n black.
  n_pairs <- nrow(ped1)
  counts <- numeric(n_pairs)
  remaining_offspring <- total_offspring
  remaining_sample <- num_offspring

  for (i in seq_len(n_pairs)) {
    if (remaining_sample <= 0) break
    pair_off <- ped1$off[i]
    remaining_after <- remaining_offspring - pair_off
    if (remaining_after <= 0) {
      # Last pair holding any offspring: it absorbs the remaining sample
      counts[i] <- remaining_sample
    } else {
      counts[i] <- stats::rhyper(nn = 1, m = pair_off, n = remaining_after,
                                 k = remaining_sample)
    }
    remaining_offspring <- remaining_after
    remaining_sample <- remaining_sample - counts[i]
  }

  ped1$off1 <- counts
  ped1$dads <- as.character(ped1$dads)
  ped1$moms <- as.character(ped1$moms)

  # Order rows by sampling weight and fix the column order
  ped1 <- ped1[order(ped1$probs), c("mp", "dads", "moms", "off", "probs", "off1")]

  return(ped1)
}
