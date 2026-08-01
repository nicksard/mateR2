#' Extract High-Fidelity Posterior Ensemble of Mate-Pair Count Vectors
#'
#' @param mcmc_res List object returned by generate_map_table().
#' @param target_Np Numeric. Your target for total parents.
#' @param target_SR Numeric. Your target sex ratio.
#' @param target_MM Numeric. Your target mean mates.
#' @param burn_in Integer. The burn-in used during the MCMC run.
#' @param thin Integer. The thinning interval used during the MCMC run.
#' @param n_samples Integer. Number of posterior draws to retain (default: 100).
#' @param max_error_pct Numeric. Max allowable relative Euclidean error (default: 0.05).
#' @export
sample_posterior_ensemble <- function(mcmc_res, target_Np, target_SR, target_MM,
                                      burn_in, thin,
                                      n_samples = 100, max_error_pct = 0.05) {

  # 1. Safely locate history and samples in the nested list structure
  history <- if (!is.null(mcmc_res$mcmc_output$history)) mcmc_res$mcmc_output$history else mcmc_res$history
  samples <- if (!is.null(mcmc_res$samples)) mcmc_res$samples else mcmc_res$mcmc_output$samples

  if (is.null(history) || is.null(samples)) {
    stop("Could not locate 'history' or 'samples' inside the provided MCMC object.")
  }

  # 2. Reconstruct the exact C++ saving indices to align history with the thinned samples
  # The C++ saves when: t >= burn_in && ((t + 1 - burn_in) % thin == 0)
  n_iter <- max(history$iteration, na.rm = TRUE)
  saved_indices <- seq(burn_in + thin, n_iter, by = thin)

  # Safety check for dimension alignment
  if (length(saved_indices) != length(samples)) {
    min_len <- min(length(saved_indices), length(samples))
    saved_indices <- saved_indices[1:min_len]
    samples <- samples[1:min_len]
  }

  # Subset history to ONLY the rows that correspond to our saved count vectors
  matched_history <- history[saved_indices, ]

  # 3. Calculate Relative Euclidean Error Distance
  rel_error <- sqrt(
    ((matched_history$Np - target_Np) / target_Np)^2 +
      ((matched_history$sr - target_SR) / target_SR)^2 +
      ((matched_history$mean_mates - target_MM) / target_MM)^2
  )

  matched_history$target_dist <- rel_error

  # 4. Filter for high-fidelity draws
  valid_rows <- which(rel_error <= max_error_pct)

  if (length(valid_rows) < n_samples) {
    warning(sprintf(
      "Only %d samples met the %.1f%% error cutoff. Retaining the top %d samples with lowest error.",
      length(valid_rows), max_error_pct * 100, min(n_samples, length(rel_error))
    ))
    # Grab indices with the lowest absolute error
    selected_valid_rows <- order(rel_error)[1:min(n_samples, length(rel_error))]
  } else {
    # Thin evenly across the passing subset to reduce sequential autocorrelation
    selected_valid_rows <- valid_rows[round(seq(1, length(valid_rows), length.out = n_samples))]
  }

  # 5. Extract the MAP state
  if (!is.null(mcmc_res$map_table)) {
    map_c_k <- mcmc_res$map_table
  } else {
    # Fallback if map_table isn't explicitly saved
    best_idx <- which.max(matched_history$log_prob)
    map_c_k <- samples[[best_idx]]
  }

  # 6. Bundle Output
  res <- list(
    map_c_k      = map_c_k,
    ensemble_c_k = samples[selected_valid_rows],
    diagnostics  = matched_history[selected_valid_rows, ],
    targets      = list(Np = target_Np, SR = target_SR, MM = target_MM)
  )

  class(res) <- "mateR2_ensemble"
  return(res)
}
