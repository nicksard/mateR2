#' @title Run and Analyze the MCMC Breeding Matrix Sampler
#' @description This is the main wrapper function that orchestrates the entire MCMC
#'   simulation. It sets up the initial conditions, runs the core C++ sampler,
#'   analyzes the output to find the Maximum a Posteriori (MAP) estimate, and
#'   calculates key statistics for the resulting breeding matrix.
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
#' @param initial_method A character string "simple" to use the default 1:1
#'   initialization, or a numeric vector representing a custom starting state.
#' @param seed An optional seed for the random number generator to ensure
#'   reproducibility.
#' @return A list containing the MAP estimate table, the summary statistics
#'   for the MAP estimate, and the raw MCMC output.
#' @import Rcpp
#' @import RcppProgress
#' @export
generate_map_table <- function(
    Np_target, sr_target, mean_mates_target,
    max_males_per_female, max_females_per_male,
    decay_constant = -0.5, np_weight = 10.0, sr_weight = 1.0, mm_weight = 10.0,
    n_iter = 200000, burn_in = 20000, thin = 20,
    initial_method = "simple", seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  print("--- Setting up MCMC ---")
  config_info <- create_config_info(max_males_per_female, max_females_per_male)
  
  if (is.character(initial_method) && initial_method == "simple") {
    initial_counts <- create_initial_counts(config_info, Np_target, sr_target)
  } else if (is.numeric(initial_method) && length(initial_method) == nrow(config_info)){
    initial_counts <- initial_method
    print("Using provided initial_counts vector.")
  } else { stop("Invalid initial_method provided.") }
  
  target_values <- list(Np_target = Np_target, sr_target = sr_target, mean_mates_target = mean_mates_target)
  mcmc_params <- list(n_iter = n_iter, burn_in = burn_in, thin = thin,
                      decay_constant = decay_constant, np_weight = np_weight,
                      sr_weight = sr_weight, mm_weight = mm_weight)
  
  print("--- Running MCMC Sampler ---")
  run_time <- system.time({
    output <- run_mcmc_sampler_cpp(
      initial_counts = initial_counts, config_info = config_info,
      target_values = target_values, mcmc_params = mcmc_params
    )
  })
  print("Run Time:"); print(run_time)
  
  if (!is.list(output) || is.null(output$samples) || is.null(output$history) || is.null(output$acceptance_rate)) {
    warning("MCMC output object is not valid. Cannot find MAP estimate.")
    return(list(map_table=NULL, map_stats=NULL, mcmc_output=output))
  }
  
  print("--- Finding MAP Estimate ---")
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
    
    for (i in 1:num_samples) {
      log_prob_result <- tryCatch({
        Counts <- output$samples[[i]]
        if(is.null(Counts) || length(Counts) != length(males_vec)) stop("Invalid Counts vector in samples")
        
        Nm <- sum(males_vec * Counts); Nf <- sum(females_vec * Counts)
        Np_real <- Nm + Nf; sr_real <- ifelse(Nf > 0, Nm / Nf, Inf)
        Total_Matings <- sum(males_vec * females_vec * Counts)
        Overall_Mean_Mates_actual <- ifelse(Np_real > 0, Total_Matings / Np_real, 0)
        
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
      
      final_Nm <- sum(map_table$Males * map_table$MAP_Count)
      final_Nf <- sum(map_table$Females * map_table$MAP_Count)
      final_Np <- final_Nm + final_Nf
      final_SR <- ifelse(final_Nf > 0, final_Nm / final_Nf, Inf)
      final_TM <- sum(map_table$Males * map_table$Females * map_table$MAP_Count)
      final_MM <- ifelse(final_Np > 0, final_TM / final_Np, 0)
      
      map_stats <- list(Np = final_Np, SR = final_SR, MeanMates = final_MM,
                        MaxLogProb = best_log_prob)
      
      print("--- MAP Estimate Found ---")
      print(map_table[map_table$MAP_Count > 0, c("Block","Males","Females","MAP_Count")])
      print("Stats for MAP table:"); print(map_stats)
      
    } else {
      warning("Could not find any valid sample with finite log probability for MAP.")
    }
  } else {
    print("No samples available to find MAP estimate.")
  }
  
  print("--- Returning Output ---")
  return(list(
    map_table = map_table,
    map_stats = map_stats,
    mcmc_output = output
  ))
}
