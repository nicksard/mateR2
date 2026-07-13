#' @title Randomize a Breeding Matrix via Bipartite Edge-Swapping
#'
#' @description Implements a joint-degree-preserving bipartite randomization routine
#'   (Curveball algorithm) to rewire the connections in a binary breeding matrix.
#'   This function controls the topological entropy of the network while
#'   rigorously preserving the exact number of mates for every individual
#'   (the demographic marginal totals).
#'
#' @param binary_IIM A binary Individual-by-Individual Matrix (IIM) where '1's
#'   indicate a mating event. This matrix should be the output of
#'   `mp_table_to_matrix()`.
#' @param I A numeric mixing intensity parameter bounded between 0 and 1.
#'   At `I = 0`, the matrix retains its highly modular, isolated initial state.
#'   At `I = 1`, the matrix reaches its maximum-entropy ceiling (panmixia).
#'   Default is 1.0.
#' @param C A fixed internal scaling constant defining the panmictic
#'   upper bound of edge swaps relative to the mixing time E*ln(E). Default is 10.
#'
#' @return A new, binary IIM that is a randomized version of the input, with
#'   identical row and column sums, and structurally randomized spatial indexing.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: a simple 4x4 block matrix
#' example_mat <- matrix(c(1,1,0,0, 1,1,0,0, 0,0,1,1, 0,0,1,1), nrow = 4, ncol = 4)
#'
#' # Maximum randomization (Uniform microcanonical configuration)
#' randomized_mat <- randomize_mating_structure(example_mat, I = 1.0)
#'
#' # Minimum randomization (Preserves block-diagonals, shuffles IDs only)
#' baseline_mat <- randomize_mating_structure(example_mat, I = 0.0)
#' }
randomize_mating_structure <- function(binary_IIM, I = 1.0, C = 10) {

  # --- 1. Input Validation ---
  if (I < 0 || I > 1) {
    stop("Mixing intensity parameter 'I' must be strictly bounded between 0 and 1.")
  }

  # --- 2. Initialization & Entropy Scaling ---
  # E: Total number of active links (edges)
  E <- sum(binary_IIM)

  # Early exit for empty or trivial matrices
  if (E < 2) {
    stop("The input matrix has fewer than 2 matings (edges). Edge-swapping is impossible.")
  }

  # T_max: The global panmictic upper bound.
  # Mixing time for edge-swapping Markov chains scales as O(E * ln(E)).
  T_max <- C * E * log(E)

  # T_swaps: Absolute number of active swaps to execute based on intensity 'I'
  T_swaps <- round(I * T_max)

  S_success <- 0 # Successful swaps counter

  # Create a list of all existing edges (coordinates of 1s)
  edge_list <- which(binary_IIM == 1, arr.ind = TRUE)

  # Use a safety break to prevent infinite loops in highly dense/sparse matrices
  max_attempts <- max(10 * T_swaps, 1)
  attempt_count <- 0

  # --- 3. The Iterative Curveball Loop ---
  # If I = 0, this loop is bypassed completely, preserving the initial modular blocks
  if (T_swaps > 0) {
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

      # C. Validation Check 2 (The Checkerboard / No redundant edges)
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

    if (attempt_count >= max_attempts && T_swaps > 0) {
      warning(paste("Randomization reached maximum attempts before completing all", T_swaps, "swaps."))
    }
  }

  # --- 4. Index Shuffling (Preventing Structural ID Bias) ---
  random_row_idx <- sample(1:nrow(binary_IIM))
  random_col_idx <- sample(1:ncol(binary_IIM))

  shuffled_IIM <- binary_IIM[random_row_idx, random_col_idx, drop = FALSE]

  if(!is.null(rownames(binary_IIM))) rownames(shuffled_IIM) <- rownames(binary_IIM)[random_row_idx]
  if(!is.null(colnames(binary_IIM))) colnames(shuffled_IIM) <- colnames(binary_IIM)[random_col_idx]

  return(shuffled_IIM)
}
