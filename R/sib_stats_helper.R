#' @title Calculate Sibling Pair Statistics from a Breeding Matrix
#' @description Analyzes a realized breeding matrix to calculate the exact
#'   number of full-sibling (FS), maternal half-sibling (MHS), and paternal
#'   half-sibling (PHS) pairwise relationships. Also returns maximum family sizes.
#'   This is critical for evaluating the joint-likelihood configuration space
#'   used by pedigree reconstruction software like COLONY.
#' @param mat The breeding matrix to be analyzed. A numeric matrix where rows
#'   represent males, columns represent females, and cells represent fecundity
#'   (number of offspring).
#' @return A data frame containing the total counts of sibling pairs
#'   (`FS_pairs`, `MHS_pairs`, `PHS_pairs`) and the maximum group sizes
#'   (`max_FS_size`, `max_maternal_family`, `max_paternal_family`).
#' @examples
#' # Assume 'my_matrix' is a breeding matrix (rows=males, cols=females)
#' # my_matrix <- matrix(c(2, 0, 1, 0, 4, 3), nrow=2, byrow=TRUE)
#' # sib_statistics <- sib.stats(my_matrix)
#' # print(sib_statistics)
#' @export
sib.stats <- function(mat) {

  # --- Input Validation ---
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("Input 'mat' must be a matrix or data frame.")
  }
  mat <- as.matrix(mat)

  # --- 1. Full-Siblings (FS) ---
  # Pairs of offspring from the exact same male-female combination
  # Formula: n * (n - 1) / 2 for every cell
  fs_pairs <- sum((mat * (mat - 1)) / 2)

  # --- 2. Maternal Half-Siblings (MHS) ---
  # Offspring sharing a mother, but having different fathers
  fem_rs <- colSums(mat) # Total offspring per female
  total_maternal_pairs <- sum((fem_rs * (fem_rs - 1)) / 2)
  mhs_pairs <- total_maternal_pairs - fs_pairs

  # --- 3. Paternal Half-Siblings (PHS) ---
  # Offspring sharing a father, but having different mothers
  male_rs <- rowSums(mat) # Total offspring per male
  total_paternal_pairs <- sum((male_rs * (male_rs - 1)) / 2)
  phs_pairs <- total_paternal_pairs - fs_pairs

  # --- 4. Compile Statistics ---
  stats_df <- data.frame(
    FS_pairs  = fs_pairs,
    MHS_pairs = mhs_pairs,
    PHS_pairs = phs_pairs,

    # Maximum group sizes (Useful for understanding the "mega-families")
    max_FS_size              = max(mat),
    max_maternal_family_size = max(fem_rs),
    max_paternal_family_size = max(male_rs)
  )

  return(stats_df)
}
