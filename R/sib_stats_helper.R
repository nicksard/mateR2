#' @title Calculate Sibling Pair Statistics and Family Size Distributions
#'
#' @description Analyzes a realized breeding matrix to calculate total counts of
#'   full-sibling (FS), maternal half-sibling (MHS), and paternal half-sibling (PHS)
#'   dyads. Additionally computes offspring-level distribution metrics (mean,
#'   minimum, maximum, variance, and standard deviation) for full-sibling,
#'   maternal half-sibling, and paternal half-sibling group sizes.
#'
#' @details
#' Sibling counts are calculated at both the dyadic level (total pairs) and individual
#' juvenile level. The distribution metrics (`mean_FS`, `mean_MHS`, `var_MHS`, etc.)
#' evaluate the sibling environment from the perspective of each sampled offspring in
#' the total cohort.
#'
#' The per-offspring means are consistent with the dyad counts: `mean_FS` equals
#' `2 * FS_pairs / sum(mat)`, and likewise for the half-sibling classes. Reporting
#' the per-offspring form directly avoids having to derive it downstream.
#'
#' @param mat A numeric matrix or data frame representing the realized breeding matrix,
#'   where rows represent males (fathers), columns represent females (mothers), and cell
#'   values indicate fecundity (number of surviving offspring produced by that pair).
#'
#' @return A single-row \code{data.frame} containing:
#'   \item{FS_pairs}{Total number of full-sibling dyads.}
#'   \item{MHS_pairs}{Total number of maternal half-sibling dyads.}
#'   \item{PHS_pairs}{Total number of paternal half-sibling dyads.}
#'   \item{max_FS_size}{Maximum fecundity produced by a single mating pair.}
#'   \item{max_maternal_family_size}{Maximum total offspring produced by a female.}
#'   \item{max_paternal_family_size}{Maximum total offspring produced by a male.}
#'   \item{mean_FS}{Mean number of full-siblings per offspring.}
#'   \item{min_FS}{Minimum number of full-siblings per offspring.}
#'   \item{max_FS}{Maximum number of full-siblings per offspring.}
#'   \item{var_FS}{Variance in full-siblings per offspring.}
#'   \item{sd_FS}{Standard deviation in full-siblings per offspring.}
#'   \item{mean_MHS}{Mean number of maternal half-siblings per offspring.}
#'   \item{min_MHS}{Minimum number of maternal half-siblings per offspring.}
#'   \item{max_MHS}{Maximum number of maternal half-siblings per offspring.}
#'   \item{var_MHS}{Variance in maternal half-siblings per offspring.}
#'   \item{sd_MHS}{Standard deviation in maternal half-siblings per offspring.}
#'   \item{mean_PHS}{Mean number of paternal half-siblings per offspring.}
#'   \item{min_PHS}{Minimum number of paternal half-siblings per offspring.}
#'   \item{max_PHS}{Maximum number of paternal half-siblings per offspring.}
#'   \item{var_PHS}{Variance in paternal half-siblings per offspring.}
#'   \item{sd_PHS}{Standard deviation in paternal half-siblings per offspring.}
#'
#' @examples
#' # Define a 2 male x 3 female breeding matrix
#' my_matrix <- matrix(c(2, 0, 1, 0, 4, 3), nrow = 2, byrow = TRUE)
#'
#' # Calculate sibling pair statistics and family size distributions
#' sib_statistics <- sib.stats(my_matrix)
#' print(sib_statistics)
#'
#' @export
sib.stats <- function(mat) {

  # --- Input Validation ---
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("Input 'mat' must be a matrix or data frame.")
  }
  mat <- as.matrix(mat)

  if (sum(mat) == 0) {
    stop("Input 'mat' contains zero total offspring.")
  }

  # --- Marginal Totals ---
  male_rs <- rowSums(mat) # Total offspring per male
  fem_rs  <- colSums(mat) # Total offspring per female

  # --- 1. Total Sibling Dyads ---
  fs_pairs  <- sum((mat * (mat - 1)) / 2)
  total_maternal_pairs <- sum((fem_rs * (fem_rs - 1)) / 2)
  mhs_pairs <- total_maternal_pairs - fs_pairs

  total_paternal_pairs <- sum((male_rs * (male_rs - 1)) / 2)
  phs_pairs <- total_paternal_pairs - fs_pairs

  # --- 2. Offspring-Level Sibling Distribution Vectors ---
  # Cell-by-cell sibling counts experienced by an offspring in cell (i, j). An
  # offspring's full sibs are the others in its own cell; its maternal (paternal)
  # half sibs are its mother's (father's) other offspring, excluding those.
  fs_matrix  <- mat - 1
  mhs_matrix <- outer(rep(1, nrow(mat)), fem_rs) - mat
  phs_matrix <- outer(male_rs, rep(1, ncol(mat))) - mat

  # Expand cell counts into an offspring-level vector (weighted by cell fecundity)
  fs_vector  <- rep(as.vector(fs_matrix),  times = as.vector(mat))
  mhs_vector <- rep(as.vector(mhs_matrix), times = as.vector(mat))
  phs_vector <- rep(as.vector(phs_matrix), times = as.vector(mat))

  # --- 3. Compile Summary Statistics ---
  stats_df <- data.frame(
    # Total Dyad Counts
    FS_pairs  = fs_pairs,
    MHS_pairs = mhs_pairs,
    PHS_pairs = phs_pairs,

    # Maximum Family Group Sizes
    max_FS_size              = max(mat),
    max_maternal_family_size = max(fem_rs),
    max_paternal_family_size = max(male_rs),

    # Full-Sibling Distributions (Per Offspring)
    mean_FS = mean(fs_vector),
    min_FS  = min(fs_vector),
    max_FS  = max(fs_vector),
    var_FS  = stats::var(fs_vector),
    sd_FS   = stats::sd(fs_vector),

    # Maternal Half-Sibling Distributions (Per Offspring)
    mean_MHS = mean(mhs_vector),
    min_MHS  = min(mhs_vector),
    max_MHS  = max(mhs_vector),
    var_MHS  = stats::var(mhs_vector),
    sd_MHS   = stats::sd(mhs_vector),

    # Paternal Half-Sibling Distributions (Per Offspring)
    mean_PHS = mean(phs_vector),
    min_PHS  = min(phs_vector),
    max_PHS  = max(phs_vector),
    var_PHS  = stats::var(phs_vector),
    sd_PHS   = stats::sd(phs_vector)
  )

  return(stats_df)
}
