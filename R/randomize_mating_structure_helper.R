#' @title Randomize a Breeding Matrix While Preserving Mating Degrees
#'
#' @description This function implements the edge-swapping (or switching) algorithm
#'   to randomize the connections in a binary breeding matrix. It rigorously
#'   preserves the number of mates for every individual (the degree sequence),
#'   while creating a more fluid, panmictic population structure.
#'
#' @param binary_IIM A binary Individual-by-Individual Matrix (IIM) where '1's
#'   indicate a mating event. This matrix should be the output of
#'   `mp_table_to_matrix()`.
#' @param SPE_target The target number of swaps per edge (SPE) to perform. A
#'   higher value (e.g., 100) will lead to a more thoroughly randomized matrix.
#'
#' @return A new, binary IIM that is a randomized version of the input, with
#'   identical row and column sums (i.e., preserved degrees).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: a simple 2x2 block matrix
#' example_mat <- matrix(c(1,1,0,0, 1,1,0,0, 0,0,1,1, 0,0,1,1), nrow = 4, ncol = 4)
#' randomized_mat <- randomize_mating_structure(example_mat, SPE_target = 100)
#' }
randomize_mating_structure <- function(binary_IIM, SPE_target) {
  # --- 1. Initialization ---
  E <- sum(binary_IIM) # Total number of edges
  T_swaps <- SPE_target * E
  S_success <- 0 # Successful swaps counter

  # Create a list of all existing edges (coordinates of 1s)
  edge_list <- which(binary_IIM == 1, arr.ind = TRUE)

  # Use a safety break to prevent infinite loops
  max_attempts <- 10 * T_swaps
  attempt_count <- 0

  # --- 2. The Iterative Loop ---
  while (S_success < T_swaps && attempt_count < max_attempts) {
    attempt_count <- attempt_count + 1

    # A. Random Selection: Choose two distinct edge indices
    idx <- sample(1:nrow(edge_list), 2)
    M1 <- edge_list[idx[1], 1]
    FA <- edge_list[idx[1], 2]
    M2 <- edge_list[idx[2], 1]
    FB <- edge_list[idx[2], 2]

    # B. Validation Check 1 (Distinct Individuals)
    if (M1 == M2 || FA == FB) {
      next
    }

    # C. Validation Check 2 (The Checkerboard)
    if (binary_IIM[M1, FB] == 0 && binary_IIM[M2, FA] == 0) {

      # D. Execution: Perform the swap in the matrix
      binary_IIM[M1, FA] <- 0
      binary_IIM[M2, FB] <- 0
      binary_IIM[M1, FB] <- 1
      binary_IIM[M2, FA] <- 1

      # E. CRITICAL FIX: Update the edge list to reflect the new reality
      edge_list[idx[1], 2] <- FB
      edge_list[idx[2], 2] <- FA

      S_success <- S_success + 1
    }
  }

  if (attempt_count >= max_attempts) {
    warning("Randomization reached maximum attempts before completing all swaps.")
  }

  # --- 3. Index Shuffling (Preventing Structural ID Bias) ---
  # The initial block-diagonal construction groups individuals by block type.
  # Shuffling the rows and columns ensures that an individual's ID (index)
  # is not correlated with their mating degree (number of mates).
  random_row_idx <- sample(1:nrow(binary_IIM))
  random_col_idx <- sample(1:ncol(binary_IIM))

  shuffled_IIM <- binary_IIM[random_row_idx, random_col_idx, drop = FALSE]

  return(shuffled_IIM)
}
