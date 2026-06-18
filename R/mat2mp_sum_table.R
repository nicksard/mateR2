#' @title Convert a Breeding Matrix to a Mate Pair Summary Table
#' @description This function takes a breeding matrix and converts it into a summary
#'   data frame of mate pair counts, similar to the format used by the MCMC sampler.
#' @param mat A breeding matrix where rows are males and columns are females.
#'   The function assumes a structured or block-like matrix.
#' @return A data frame with three columns: `Males`, `Females`, and `Count`,
#'   summarizing the number of males, females, and the number of mating blocks
#'   of that type.
#' @importFrom dplyr mutate filter group_by summarize
#' @importFrom tidyr pivot_longer
#' @export
mat2mp_sum_table <- function(mat) {
  # --- Input Validation ---
  if (!is.matrix(mat) || any(mat < 0) || any(mat != 0 & mat != 1)) {
    stop("Input 'mat' must be a binary matrix (0s and 1s) with non-negative values.")
  }
  
  # --- Calculate Row and Column Sums ---
  # These represent the number of mates for each male and female
  row_sums <- rowSums(mat)
  col_sums <- colSums(mat)
  
  # --- Create a Mate Pair Table ---
  # Combine row and column sums into a data frame
  mate_pairs <- data.frame(
    Males = row_sums,
    Females = col_sums
  )
  
  # --- Summarize Counts ---
  # Group by the male and female mate counts to get a summary table
  summary_table <- mate_pairs %>%
    group_by(Males, Females) %>%
    summarize(Count = n(), .groups = "drop") %>%
    as.data.frame()
  
  return(summary_table)
}
