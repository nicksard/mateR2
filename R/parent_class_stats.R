#' @title Calculate Parent Abundance across Singleton, Doubleton, Tripleton, and Quadrupleton Classes
#'
#' @description Evaluates a realized breeding matrix to count the number of detected
#'   maternal, paternal, and overall parents that produced exactly 1 (singleton),
#'   2 (doubleton), 3 (tripleton), or 4 (quadrupleton) offspring.
#'
#' @details
#' Individual parent reproductive success (RS) is derived from row sums (paternal)
#' and column sums (maternal). Detected parents are defined as those with $RS > 0$.
#' The function returns counts of parents falling into discrete offspring yield
#' tiers, which is particularly useful for evaluating coverage and Chao/Jackknife
#' richness estimation models in pedigree accumulation studies.
#'
#' @param mat A numeric matrix or data frame representing a realized breeding matrix,
#'   where rows represent males (fathers), columns represent females (mothers), and cell
#'   values indicate fecundity (offspring produced).
#'
#' @return A \code{data.frame} with 3 rows (\code{maternal}, \code{paternal}, \code{overall})
#'   and 5 numeric columns:
#'   \item{detected}{Total number of active parents with at least 1 offspring ($RS > 0$).}
#'   \item{singletons}{Number of parents with exactly 1 offspring ($RS = 1$).}
#'   \item{doubletons}{Number of parents with exactly 2 offspring ($RS = 2$).}
#'   \item{tripletons}{Number of parents with exactly 3 offspring ($RS = 3$).}
#'   \item{quadrupletons}{Number of parents with exactly 4 offspring ($RS = 4$).}
#'
#' @examples
#' # Define example 2 male x 3 female breeding matrix
#' my_matrix <- matrix(c(2, 0, 1, 0, 4, 3), nrow = 2, byrow = TRUE)
#'
#' # Calculate parental class distributions
#' parent.class.stats(my_matrix)
#'
#' @export
parent.class.stats <- function(mat) {

  # --- Input Validation ---
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("Input 'mat' must be a matrix or data frame.")
  }
  mat <- as.matrix(mat)

  # --- Calculate Individual Reproductive Success (RS) ---
  male_rs <- rowSums(mat)  # Paternal RS
  fem_rs  <- colSums(mat)  # Maternal RS

  # Helper function to compute class counts
  summarize_rs <- function(rs_vec) {
    active_rs <- rs_vec[rs_vec > 0]
    data.frame(
      detected      = length(active_rs),
      singletons    = sum(active_rs == 1),
      doubletons    = sum(active_rs == 2),
      tripletons    = sum(active_rs == 3),
      quadrupletons = sum(active_rs == 4)
    )
  }

  # --- Compile Summary Rows ---
  maternal_df <- summarize_rs(fem_rs)
  paternal_df <- summarize_rs(male_rs)
  overall_df  <- summarize_rs(c(fem_rs, male_rs))

  # Combine into a single data frame
  res_df <- rbind(
    maternal = maternal_df,
    paternal = paternal_df,
    overall  = overall_df
  )

  return(res_df)
}
