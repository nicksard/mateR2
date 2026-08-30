#' @title Run and Analyze the MCMC Breeding Matrix Sampler
#' @description This is the main wrapper function that orchestrates the entire MCMC
#'   simulation. It sets up the initial conditions, runs the core C++ sampler,
#'   analyzes the output to find the Maximum a Posteriori (MAP) estimate, and
#'   calculates key statistics for the resulting breeding matrix. Optionally extracts
#'   a high-fidelity posterior ensemble.
#' @param Np_target The target total number of parents (males + females).
#' @param sr_target The target sex ratio (males / females).
#' @param mean_mates_target The target mean number of mates per individual.
#' @param max_males_per_female The maximum number of males a female can mate with.
#' @param max_females_per_male The maximum number of females a male can mate with.
#' @param decay_constant A negative value that penalizes complex mating structures.
#'   A larger negative value results in a higher penalty.
#' @param np_weight The weight given to the target Number of Parents (Np).
#' @param sr_weight The weight given to the target Sex Ratio (SR).
#' @param mm_weight The weight given to the target Mean Mates (MM).
#' @param n_iter The total number of iterations for the MCMC sampler.
#' @param burn_in The number of initial iterations to discard before sampling begins.
#' @param thin The thinning interval for collecting samples from the MCMC chain.
#' @param initial_method A character string ("auto") to dynamically calculate the
#'   optimal starting state based on targets, or a numeric vector representing a
#'   custom starting state of block counts.
#' @param seed An optional seed for the random number generator to ensure reproducibility.
#' @param sample_ensemble Logical. If TRUE, extracts a thinned posterior ensemble matching the demographic targets.
#' @param n_ensemble Integer. The number of posterior draws to retain in the ensemble (default: 100).
#' @param max_error_pct Numeric. The maximum allowable relative Euclidean distance from targets (default: 0.05).
#' @param verbose Logical. If TRUE, prints progress and MAP outputs to the console (default: FALSE).
#' @param show_progress Logical. If TRUE, the C++ sampler draws a progress bar
#'   and prints the final acceptance rate. Defaults to \code{interactive()}, so
#'   scripted and parallel runs are silent; a progress bar written from a
#'   parallel worker is noise in a log and can stall console I/O on Windows.
#' @return A list containing the MAP estimate table, the summary statistics for the MAP estimate,
#'   the raw MCMC output, and optionally the high-fidelity posterior ensemble.
#' @import Rcpp
#' @import RcppProgress
#' @export
generate_map_table <- function(
    Np_target, sr_target, mean_mates_target,
    max_males_per_female, max_females_per_male,
    decay_constant = -0.5, np_weight = 10.0, sr_weight = 1.0, mm_weight = 10.0,
    n_iter = 200000, burn_in = 20000, thin = 20,
    initial_method = "auto", seed = NULL,
    sample_ensemble = FALSE, n_ensemble = 100, max_error_pct = Inf, verbose = FALSE,
    show_progress = interactive(), complexity = c("sum", "cyclomatic", "quadratic", "asymmetry")
) {
  complexity <- match.arg(complexity)
  if (!is.null(seed)) set.seed(seed)

  if (verbose) print("--- Validating Demographic Targets ---")
  check_target_viability(sr_target, mean_mates_target, max_males_per_female, max_females_per_male)

  if (verbose) print("--- Setting up MCMC ---")
  config_info <- create_config_info(max_males_per_female, max_females_per_male,
                                    complexity = complexity)

  # THE "AUTO" VS "CUSTOM" LOGIC
  if (is.character(initial_method) && initial_method == "auto") {
    initial_counts <- create_initial_counts(config_info, Np_target, sr_target)
    if (verbose) print("Using automatically generated 'Warm Start' based on target Sex Ratio.")
  } else if (is.numeric(initial_method) && length(initial_method) == nrow(config_info)){
    initial_counts <- initial_method
    if (verbose) print("Using user-provided custom initial_counts vector.")
  } else { stop("Invalid initial_method provided.") }

  target_values <- list(Np_target = Np_target, sr_target = sr_target, mean_mates_target = mean_mates_target)
  mcmc_params <- list(n_iter = n_iter, burn_in = burn_in, thin = thin,
                      decay_constant = decay_constant, np_weight = np_weight,
                      sr_weight = sr_weight, mm_weight = mm_weight,
                      show_progress = isTRUE(show_progress))

  if (verbose) print("--- Running MCMC Sampler ---")
  run_time <- system.time({
    output <- run_mcmc_sampler_cpp(
      initial_counts = initial_counts, config_info = config_info,
      target_values = target_values, mcmc_params = mcmc_params
    )
  })
  if (verbose) { print("Run Time:"); print(run_time) }

  if (!is.list(output) || is.null(output$samples) || is.null(output$history) || is.null(output$acceptance_rate)) {
    warning("MCMC output object is not valid. Cannot find MAP estimate.")
    return(list(map_table=NULL, map_stats=NULL, mcmc_output=output))
  }

  if (verbose) print("--- Processing Samples (MAP & Ensemble) ---")
  map_sample_counts <- NULL
  map_stats <- list(Np=NA, SR=NA, MeanMates=NA, MaxLogProb=NA)
  map_table <- config_info
  map_table$MAP_Count <- NA_real_

  num_samples <- length(output$samples)

  if (num_samples > 0) {
    map_log_probs <- numeric(num_samples)
    males_vec <- config_info$Males
    females_vec <- config_info$Females
    comp_diff_vec <- config_info$Complexity_Diff

    # --- Data Tracking for Ensemble Extraction ---
    ens_Np <- numeric(num_samples)
    ens_SR <- numeric(num_samples)
    ens_MM <- numeric(num_samples)

    for (i in 1:num_samples) {
      log_prob_result <- tryCatch({
        Counts <- output$samples[[i]]
        if(is.null(Counts) || length(Counts) != length(males_vec)) stop("Invalid Counts vector in samples")

        Nm <- sum(males_vec * Counts); Nf <- sum(females_vec * Counts)
        Np_real <- Nm + Nf; sr_real <- ifelse(Nf > 0, Nm / Nf, Inf)

        Total_Matings <- sum(males_vec * females_vec * Counts)

        # --- FIXED BIOLOGICAL MEAN MATES ALGEBRA ---
        male_mean_mates <- ifelse(Nm > 0, Total_Matings / Nm, 0)
        female_mean_mates <- ifelse(Nf > 0, Total_Matings / Nf, 0)
        Overall_Mean_Mates_actual <- (male_mean_mates + female_mean_mates) / 2

        # Track statistics dynamically for the ensemble filtering
        ens_Np[i] <- Np_real
        ens_SR[i] <- sr_real
        ens_MM[i] <- Overall_Mean_Mates_actual

        Np_score <- calculate_closeness_score(Np_real, Np_target)
        sr_score <- calculate_closeness_score(sr_real, sr_target)
        mm_score <- calculate_closeness_score(Overall_Mean_Mates_actual, mean_mates_target)

        current_log_prob <- -Inf
        if (Np_score >= .Machine$double.eps && sr_score >= .Machine$double.eps && mm_score >= .Machine$double.eps) {
          log_score_part <- np_weight * log(Np_score) + sr_weight * log(sr_score) + mm_weight * log(mm_score)
          total_complexity_diff <- sum(comp_diff_vec * Counts)
          log_decay_part <- decay_constant * total_complexity_diff
          current_log_prob <- log_score_part + log_decay_part
        }
        current_log_prob
      }, error = function(e) {
        warning(paste("Error calculating log prob for sample index", i, ":", e$message))
        return(-Inf)
      })
      map_log_probs[i] <- log_prob_result
    }

    finite_log_probs <- map_log_probs[is.finite(map_log_probs)]
    if (length(finite_log_probs) > 0) {
      best_log_prob <- max(finite_log_probs)
      best_sample_index <- which(map_log_probs == best_log_prob)[1]

      map_sample_counts <- output$samples[[best_sample_index]]
      map_table$MAP_Count <- map_sample_counts

      map_stats <- list(
        Np = ens_Np[best_sample_index],
        SR = ens_SR[best_sample_index],
        MeanMates = ens_MM[best_sample_index],
        MaxLogProb = best_log_prob
      )

      if (verbose) {
        print("--- MAP Estimate Found ---")
        print(map_table[map_table$MAP_Count > 0, c("Block","Males","Females","MAP_Count")])
        print("Stats for MAP table:"); print(map_stats)
      }

    } else {
      warning("Could not find any valid sample with finite log probability for MAP.")
    }

    # =====================================================================
    # --- POSTERIOR ENSEMBLE SAMPLING LOGIC ---
    # =====================================================================
    ensemble_c_k <- NULL
    ensemble_diags <- NULL

    if (sample_ensemble) {
      if (verbose) print("--- Extracting High-Fidelity Ensemble ---")

      # Calculate Relative Euclidean Distance for all thinned samples
      rel_error <- sqrt(
        ((ens_Np - Np_target) / Np_target)^2 +
          ((ens_SR - sr_target) / sr_target)^2 +
          ((ens_MM - mean_mates_target) / mean_mates_target)^2
      )

      valid_idx <- which(rel_error <= max_error_pct)

      if (length(valid_idx) < n_ensemble) {
        if (verbose) warning(sprintf("Only %d samples met the %.1f%% error cutoff. Retaining top %d.", length(valid_idx), max_error_pct * 100, min(n_ensemble, length(rel_error))))
        selected_idx <- order(rel_error)[1:min(n_ensemble, length(rel_error))]
      } else {
        # Thin evenly across the valid subset to reduce autocorrelation
        selected_idx <- valid_idx[round(seq(1, length(valid_idx), length.out = n_ensemble))]
      }

      ensemble_c_k <- output$samples[selected_idx]
      ensemble_diags <- data.frame(
        Sample_Index = selected_idx,
        Np = ens_Np[selected_idx],
        SR = ens_SR[selected_idx],
        MeanMates = ens_MM[selected_idx],
        Target_Error = rel_error[selected_idx]
      )
    }

  } else {
    print("No samples available to find MAP estimate.")
  }

  if (verbose) print("--- Returning Output ---")

  # Final Return Object
  res <- list(
    map_table   = map_table,
    map_stats   = map_stats,
    mcmc_output = output
  )

  # Append ensemble if requested
  if (sample_ensemble) {
    res$ensemble_c_k   <- ensemble_c_k
    res$ensemble_diags <- ensemble_diags
  }

  return(res)
}
