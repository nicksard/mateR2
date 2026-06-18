#' @title Convert a Pedigree Data Frame to a Breeding Matrix
#' @description This function takes a pedigree in a long data frame format and
#'   converts it into a wide breeding matrix (Individual-by-Individual Matrix),
#'   where rows represent males (dads) and columns represent females (moms).
#' @param ped A data frame representing a pedigree, with columns `dad` and `mom`.
#' @return A matrix where the cell at `[i, j]` contains the number of offspring
#'   produced by male `i` and female `j`.
#' @importFrom stats xtabs
#' @export
ped2mat <- function(ped) {
  # Create a contingency table using xtabs to count offspring per mom-dad pair
  mat <- xtabs(~ mom + dad, data = ped)
  
  # Convert to a regular matrix and transpose
  # Transposing makes sure dads are rows and moms are columns
  mat <- as.matrix(mat)
  mat <- t(mat)
  
  # Ensure dimensions are named consistently for clarity, even if dropped later
  dimnames(mat) <- list(dads = rownames(mat), moms = colnames(mat))
  
  return(mat)
}
