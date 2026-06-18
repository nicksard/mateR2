#' @title Convert MCMC Mate Pair Table to a Breeding Matrix
#' @description This function takes a summary table of breeding configurations
#'   (e.g., the output from a MCMC run) and converts it into a structured, binary
#'   Individual-by-Individual Matrix (IIM).
#' @param mp_table A data frame with at least the columns `Males`, `Females`, and `Count`,
#'   representing the number of males, females, and the number of blocks of that
#'   type.
#' @return A binary matrix where rows represent males and columns represent females.
#'   A value of 1 indicates a mating occurred. The matrix is block-diagonal.
#' @export
mp_table_to_matrix <- function(mp_table) {
  # --- Input Validation ---
  if (!is.data.frame(mp_table) || !all(c("Males", "Females", "MAP_Count") %in% names(mp_table))) {
    stop("Input 'mp_table' must be a data frame with columns 'Males', 'Females', and 'Count'.")
  }

  # --- Calculate Matrix Dimensions ---
  num_males <- sum(mp_table$Males * mp_table$MAP_Count)
  num_females <- sum(mp_table$Females * mp_table$MAP_Count)

  if (num_males == 0 || num_females == 0) {
    warning("Input table results in a matrix with zero dimensions. Returning an empty matrix.")
    return(matrix(0, nrow = num_males, ncol = num_females))
  }

  # --- Initialize an empty matrix ---
  mat <- matrix(0, nrow = num_males, ncol = num_females)

  # --- Populate the Matrix ---
  current_row <- 1
  current_col <- 1

  # Iterate through each row of the input table
  for (i in 1:nrow(mp_table)) {
    males <- mp_table$Males[i]
    females <- mp_table$Females[i]
    count <- mp_table$MAP_Count[i]

    # Only process if the count is greater than zero
    if (count > 0) {
      # Iterate for the number of blocks of this type
      for (j in 1:count) {
        # Check if enough space is available in the matrix
        if (current_row + males - 1 > nrow(mat) || current_col + females - 1 > ncol(mat)) {
          stop(paste0("Not enough space in the matrix to accommodate the '", males, ":", females, "' mating block."))
        }

        # Fill the sub-matrix with 1s
        mat[current_row:(current_row + males - 1), current_col:(current_col + females - 1)] <- 1

        # Update the start indices for the next block
        current_row <- current_row + males
        current_col <- current_col + females
      }
    }
  }

  return(mat)
}
