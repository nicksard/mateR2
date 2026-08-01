#' Convert a Pedigree Data Frame to a Breeding Matrix
#'
#' Takes a pedigree in long data frame format and converts it into a wide
#' breeding matrix (Individual-by-Individual Matrix), where rows represent
#' males (dads) and columns represent females (moms).
#'
#' @param ped A data frame representing a pedigree, with columns for parents
#'   (`mom`/`Mom` and `dad`/`Dad`).
#' @return A matrix where cell [i, j] contains the number of offspring
#'   produced by male i and female j.
#' @importFrom stats xtabs
#' @export
ped2mat <- function(ped) {
  if (!is.data.frame(ped)) {
    stop("Input 'ped' must be a data frame.")
  }

  # Standardize column names to lowercase for robust xtabs matching
  colnames(ped) <- tolower(colnames(ped))

  if (!all(c("mom", "dad") %in% colnames(ped))) {
    stop("Pedigree data frame must contain 'mom' (or 'Mom') and 'dad' (or 'Dad') columns.")
  }

  # Create contingency table of offspring per mom-dad pair
  mat <- xtabs(~ mom + dad, data = ped)

  # Convert to regular matrix and transpose so rows = dads (males), cols = moms (females)
  mat <- t(as.matrix(mat))

  # Set dimension names
  dimnames(mat) <- list(dads = rownames(mat), moms = colnames(mat))

  return(mat)
}
