#' @title Convert a Breeding Matrix to a Pedigree Data Frame
#' @description This function takes a breeding matrix and converts it into a
#'   long-format pedigree data frame, where each row represents a unique offspring
#'   and its parents. This is the inverse of `ped2mat()`.
#' @param mat A breeding matrix where rows are males and columns are females.
#' @return A data frame with three columns: `off` (offspring ID), `mom` (mother ID),
#'   and `dad` (father ID).
#' @export
mat2ped <- function(mat) {
  # Get dimensions of the matrix
  num_rows <- nrow(mat)
  num_cols <- ncol(mat)
  
  # Check if the matrix is valid
  if (!is.matrix(mat) || !is.numeric(mat) || any(mat < 0)) {
    stop("Input 'mat' must be a numeric matrix with non-negative values.")
  }
  
  # Create a data frame with offspring, mother, and father IDs,
  # repeating rows based on the count in the matrix.
  ped <- data.frame(
    mom = rep(1:num_cols, each = num_rows),
    dad = rep(1:num_rows, times = num_cols)
  )
  
  # Repeat the rows based on the number of offspring in the matrix cells
  ped <- ped[as.vector(mat) > 0, ]
  ped_with_counts <- data.frame(
    mom = ped$mom,
    dad = ped$dad,
    offspring_count = as.vector(mat)[as.vector(mat) > 0]
  )
  
  # Expand the data frame based on offspring counts
  ped_expanded <- ped_with_counts[rep(1:nrow(ped_with_counts), ped_with_counts$offspring_count), c("mom", "dad")]
  
  # Create unique offspring IDs
  ped_expanded$off <- paste0("off_", 1:nrow(ped_expanded))
  
  # Rename columns for clarity
  ped_expanded$mom <- paste0("mom", ped_expanded$mom)
  ped_expanded$dad <- paste0("dad", ped_expanded$dad)
  
  # Reorder columns to the standard pedigree format
  ped_expanded <- ped_expanded[, c("off", "mom", "dad")]
  
  return(ped_expanded)
}
