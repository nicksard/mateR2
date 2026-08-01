#' Convert Mate-Pair Summary Table to Binary Mating Matrix
#'
#' Expands a 3-column mate-pair summary table (block male size, block female size, count)
#' into a full binary individual-by-individual matrix (IIM).
#'
#' @param mp_table Data frame containing 'Males', 'Females', and a count column
#'   ('MAP_Count', 'Count', or 'n').
#'
#' @return Binary integer matrix where rows represent individual males and columns represent females.
#' @export
mp_table_to_matrix <- function(mp_table) {
  if (!is.data.frame(mp_table)) {
    stop("Input 'mp_table' must be a data frame.")
  }

  if (!all(c("Males", "Females") %in% names(mp_table))) {
    stop("Input 'mp_table' must contain 'Males' and 'Females' columns.")
  }

  # Support flexible count column naming
  count_col <- intersect(c("MAP_Count", "Count", "n"), names(mp_table))
  if (length(count_col) == 0) {
    stop("Input 'mp_table' must contain a count column named 'MAP_Count', 'Count', or 'n'.")
  }
  count_var <- count_col[1]
  counts <- as.numeric(mp_table[[count_var]])

  num_males <- sum(mp_table$Males * counts)
  num_females <- sum(mp_table$Females * counts)

  if (num_males == 0 || num_females == 0) {
    warning("Input table results in a matrix with zero dimensions. Returning an empty matrix.")
    return(matrix(0, nrow = num_males, ncol = num_females))
  }

  mat <- matrix(0L, nrow = num_males, ncol = num_females)
  current_row <- 1
  current_col <- 1

  for (i in seq_len(nrow(mp_table))) {
    males <- mp_table$Males[i]
    females <- mp_table$Females[i]
    count <- counts[i]

    if (count > 0) {
      for (j in seq_len(count)) {
        if (current_row + males - 1 > nrow(mat) || current_col + females - 1 > ncol(mat)) {
          stop(paste0("Not enough space in the matrix to accommodate the '",
                      males, ":", females, "' mating block."))
        }
        mat[current_row:(current_row + males - 1), current_col:(current_col + females - 1)] <- 1L
        current_row <- current_row + males
        current_col <- current_col + females
      }
    }
  }

  return(mat)
}
