#' Simulate Pedigrees Across a Posterior Ensemble of Count Vectors
#'
#' Takes a list of Mate-Pair Summary Table count vectors (\eqn{\vec{c}_k}), expands each
#' into an Individual-by-Individual Matrix (IIM), applies Curveball degree-preserving
#' rewiring, populates offspring fecundities, and sub-samples juveniles to produce ready-to-use pedigrees.
#'
#' @param ensemble Either an object of class \code{"mateR2_ensemble"} returned by
#'   \code{\link{sample_posterior_ensemble}}, or a direct list of \eqn{\vec{c}_k} count vectors.
#' @param mixing_I Numeric. Curveball rewiring intensity between 0.0 (block-diagonal)
#'   and 1.0 (full panmictic mixing). Default is 0.25.
#' @param min_fecundity Integer. Minimum offspring allocated per active pair-bond (default: 20).
#' @param max_fecundity Integer. Maximum offspring allocated per active pair-bond (default: 100).
#' @param fecundity_type Character. Fitness distribution type for \code{\link{brd.mat.fitness}}
#'   (e.g., \code{"uniform"}, \code{"lognormal"}). Default is \code{"uniform"}.
#' @param juvenile_sample_size Integer. Number of juveniles to sub-sample per pedigree (default: 500).
#' @param include_map Logical. If \code{ensemble} is a \code{"mateR2_ensemble"} object,
#'   whether to include the MAP pedigree as the first element of the output list. Default is \code{TRUE}.
#'
#' @return A list of ground-truth pedigree data frames with columns \code{c("Offspring", "Mom", "Dad")}.
#'
#' @examples
#' \dontrun{
#' mcmc_res <- generate_map_table(Np_target = 500, sr_target = 1.5, mean_mates_target = 2.0)
#' ensemble <- sample_posterior_ensemble(mcmc_res, n_samples = 50)
#' pedigree_list <- simulate_pedigree_ensemble(ensemble, mixing_I = 0.25, juvenile_sample_size = 500)
#' }
#' @export
simulate_pedigree_ensemble <- function(ensemble,
                                       mixing_I = 0.25,
                                       min_fecundity = 20,
                                       max_fecundity = 100,
                                       fecundity_type = "uniform",
                                       juvenile_sample_size = 500,
                                       include_map = TRUE) {

  # --- 1. Input Parsing ---
  if (inherits(ensemble, "mateR2_ensemble")) {
    c_k_list <- ensemble$ensemble_c_k
    if (include_map && !is.null(ensemble$map_c_k)) {
      c_k_list <- c(list(MAP = ensemble$map_c_k), c_k_list)
    }
  } else if (is.list(ensemble)) {
    c_k_list <- ensemble
  } else {
    stop("Input 'ensemble' must be a list of c_k vectors or a 'mateR2_ensemble' object.")
  }

  # --- 2. Process Each Vector Through Stages 1 -> 2 -> 3 ---
  pedigrees <- lapply(seq_along(c_k_list), function(i) {
    c_k <- c_k_list[[i]]

    # Step 1: Expand c_k into block-diagonal binary IIM
    binary_mat <- mp_table_to_matrix(c_k)

    # Step 2: Stage 2 Curveball Edge Swapping
    if (mixing_I > 0) {
      mixed_mat <- randomize_mating_structure(binary_mat, intensity = mixing_I)
    } else {
      mixed_mat <- binary_mat
    }

    # Step 3: Populate Pair Offspring Fecundity
    fecund_mat <- brd.mat.fitness(
      mat      = mixed_mat,
      min.fert = min_fecundity,
      max.fert = max_fecundity,
      type     = fecundity_type
    )

    # Step 4: Sub-sample Juveniles (Field Sampling Filter)
    sampled_df <- mat.sub.sample(fecund_mat, num_offspring = juvenile_sample_size)

    # Step 5: Convert to Standard Pedigree Data Frame
    ped <- convert2ped(sampled_df)
    return(ped)
  })

  # Preserve names if present (e.g., MAP)
  if (!is.null(names(c_k_list))) {
    names(pedigrees) <- names(c_k_list)
  }

  return(pedigrees)
}
