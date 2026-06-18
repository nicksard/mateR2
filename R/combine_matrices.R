#' @title Combine Two Breeding Matrices
#' @description This function takes two separate breeding matrices and combines
#'   them into a single, larger matrix by placing them on the diagonal. This is
#'   useful for joining breeding structures from different populations or scenarios.
#' @param mat1 The first breeding matrix.
#' @param mat2 The second breeding matrix.
#' @return A new, combined matrix that contains both input matrices on the diagonal.
#' @export
combine_matrices <- function(mat1, mat2) {
  # --- Input Validation ---
  if (!is.matrix(mat1) || !is.matrix(mat2)) {
    stop("Inputs 'mat1' and 'mat2' must be matrices.")
  }
  
  # --- Create a new, empty matrix with the combined dimensions ---
  combined_mat <- matrix(0,
                         nrow = nrow(mat1) + nrow(mat2),
                         ncol = ncol(mat1) + ncol(mat2)
  )
  
  # --- Place the input matrices into the combined matrix ---
  # Place the first matrix in the top-left corner
  combined_mat[1:nrow(mat1), 1:ncol(mat1)] <- mat1
  
  # Place the second matrix in the bottom-right corner
  combined_mat[
    (nrow(mat1) + 1):(nrow(mat1) + nrow(mat2)),
    (ncol(mat1) + 1):(ncol(mat1) + ncol(mat2))
  ] <- mat2
  
  return(combined_mat)
}
