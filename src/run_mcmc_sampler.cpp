#include <Rcpp.h>
#include <vector>
#include <cmath> // For log, exp, abs, floor
#include <numeric> // For std::accumulate (optional)
#include <random> // For C++11 random number generation
#include <chrono> // For seeding

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp> // For progress bar

using namespace Rcpp;

// --- C++ Helper Function: Closeness Score ---
// @title Calculate a Closeness Score
 // @description Calculates min(Ratio, 1/Ratio), ensuring the result is > 0 for log-stability.
 // @param actual The actual value.
 // @param target_value The target value.
 // @return A double between 0 and 1 representing the closeness of the actual value to the target.
 double cpp_closeness_score(double actual, double target_value) {
   if (target_value <= 0 || !std::isfinite(target_value)) {
     return 0.0;
   }
   if (actual <= 0 || !std::isfinite(actual)) {
     return 0.0;
   }

   double ratio = actual / target_value;
   double inv_ratio = target_value / actual;
   double score = std::min(ratio, inv_ratio);

   // Return score, ensuring it's slightly above zero if calculated as zero
   return std::max(score, std::numeric_limits<double>::epsilon());
 }

// --- C++ Helper Function: Calculate Target Log Probability ---
// @title Calculate Target Log Probability
 // @description Calculates the log probability of a given breeding configuration.
 // @param Counts A numeric vector of block counts.
 // @param Males A numeric vector of males per block type.
 // @param Females A numeric vector of females per block type.
 // @param Complexity_Diff A numeric vector of complexity difference per block type.
 // @param Np_target The target number of parents.
 // @param sr_target The target sex ratio.
 // @param mean_mates_target The target mean number of mates.
 // @param decay_constant A negative value to penalize complexity.
 // @param np_weight The weight for the number of parents score.
 // @param sr_weight The weight for the sex ratio score.
 // @param mm_weight The weight for the mean mates score.
 // @return The log probability (a double).
 double cpp_calculate_target_log_prob(
     const NumericVector& Counts,
     const NumericVector& Males,
     const NumericVector& Females,
     const NumericVector& Complexity_Diff,
     double Np_target,
     double sr_target,
     double mean_mates_target,
     double decay_constant, // negative
     double np_weight,
     double sr_weight,
     double mm_weight
 ) {
   int n_types = Counts.length();
   double Nm = 0.0, Nf = 0.0, Total_Matings = 0.0, total_complexity_diff = 0.0;

   // Calculate sums using a loop
   for (int i = 0; i < n_types; ++i) {
     if (Counts[i] > 0) {
       Nm += Males[i] * Counts[i];
       Nf += Females[i] * Counts[i];
       Total_Matings += Males[i] * Females[i] * Counts[i];
       total_complexity_diff += Complexity_Diff[i] * Counts[i];
     }
   }

   double Np_real = Nm + Nf;
   double sr_real = (Nf > 0) ? (Nm / Nf) : R_PosInf;
   double Overall_Mean_Mates_actual = (Np_real > 0) ? ((2.0 * Total_Matings) / Np_real) : 0.0;

   // Calculate individual scores
   double Np_score = cpp_closeness_score(Np_real, Np_target);
   double sr_score = cpp_closeness_score(sr_real, sr_target);
   double mm_score = cpp_closeness_score(Overall_Mean_Mates_actual, mean_mates_target);

   // Check for invalid scores
   if (Np_score == 0 || sr_score == 0 || mm_score == 0) {
     return R_NegInf;
   }

   // Combine scores (log space) using weights
   double log_score_part = np_weight * std::log(Np_score) +
     sr_weight * std::log(sr_score) +
     mm_weight * std::log(mm_score);

   // Calculate Decay Part
   double log_decay_part = decay_constant * total_complexity_diff;

   return log_score_part + log_decay_part;
 }

// --- Main MCMC Sampler Function (Exported to R) ---
// @title Run MCMC Sampler
 // @description The main MCMC sampler function implemented in C++ for performance.
 // @param initial_counts A numeric vector representing the starting breeding configuration.
 // @param config_info A data frame with information on each block type.
 // @param target_values A list of target demographic parameters.
 // @param mcmc_params A list of MCMC tuning parameters.
 // @return A list containing the `samples`, `history`, and `acceptance_rate` of the MCMC run.
 // [[Rcpp::export]]
 List run_mcmc_sampler_cpp(
     NumericVector initial_counts,
     DataFrame config_info,
     List target_values,
     List mcmc_params
 ) {
   // --- Parameter Unpacking & Setup ---
   int n_iter = as<int>(mcmc_params["n_iter"]);
   int burn_in = as<int>(mcmc_params["burn_in"]);
   int thinning_interval = as<int>(mcmc_params["thin"]);
   double decay_constant = as<double>(mcmc_params["decay_constant"]);
   double np_weight = as<double>(mcmc_params["np_weight"]);
   double sr_weight = as<double>(mcmc_params["sr_weight"]);
   double mm_weight = as<double>(mcmc_params["mm_weight"]);
   double Np_target = as<double>(target_values["Np_target"]);
   double sr_target = as<double>(target_values["sr_target"]);
   double mean_mates_target = as<double>(target_values["mean_mates_target"]);

   // Extract vectors from config_info DataFrame
   NumericVector Males = config_info["Males"];
   NumericVector Females = config_info["Females"];
   NumericVector Complexity_Diff = config_info["Complexity_Diff"];
   int num_block_types = initial_counts.length();

   // Input Checks (basic examples)
   if (Males.length() != num_block_types || Females.length() != num_block_types || Complexity_Diff.length() != num_block_types) {
     stop("Config info vectors length must match initial_counts length.");
   }
   if (burn_in >= n_iter) {
     stop("burn_in must be less than n_iter.");
   }

   // --- Initialization ---
   NumericVector Current_Counts = clone(initial_counts);
   double Current_Log_Prob = cpp_calculate_target_log_prob(
     Current_Counts, Males, Females, Complexity_Diff,
     Np_target, sr_target, mean_mates_target,
     decay_constant, np_weight, sr_weight, mm_weight
   );
   if (!std::isfinite(Current_Log_Prob)) {
     stop("Initial state has -Inf log probability. Check parameters/initial state.");
   }

   // --- History & Sample Storage (Using std::vector for history) ---
   std::vector<double> history_log_prob; history_log_prob.reserve(n_iter);
   std::vector<double> history_Np; history_Np.reserve(n_iter);
   std::vector<double> history_sr; history_sr.reserve(n_iter);
   std::vector<double> history_mean_mates; history_mean_mates.reserve(n_iter);

   List results_list;
   int accepted_count = 0;

   // Progress bar setup
   Progress p(n_iter, true);

   // Use C++11 random number generation
   unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
   std::mt19937_64 gen(seed);
   std::uniform_real_distribution<double> runif_dist(0.0, 1.0);

   // --- MCMC Loop ---
   for (int t = 0; t < n_iter; ++t) {
     if (Progress::check_abort()) {
       stop("MCMC interrupted by user.");
     }

     // **a. Propose**
     NumericVector Proposed_Counts = clone(Current_Counts);
     // Sample index (0 to n-1 in C++)
     int i = std::floor(runif_dist(gen) * num_block_types);

     // Propose delta=+1 or -1
     int delta = (Proposed_Counts[i] == 0) ? 1 : ((runif_dist(gen) < 0.5) ? -1 : 1);
     Proposed_Counts[i] += delta;

     // **b. Evaluate Proposal**
     double Proposed_Log_Prob = cpp_calculate_target_log_prob(
       Proposed_Counts, Males, Females, Complexity_Diff,
       Np_target, sr_target, mean_mates_target,
       decay_constant, np_weight, sr_weight, mm_weight
     );

     // **c. Accept/Reject**
     if (std::isfinite(Proposed_Log_Prob)) {
       double log_acceptance_ratio = Proposed_Log_Prob - Current_Log_Prob;
       if (std::log(runif_dist(gen)) < log_acceptance_ratio) {
         // Accept
         Current_Counts = Proposed_Counts;
         Current_Log_Prob = Proposed_Log_Prob;
         accepted_count++;
       }
     }

     // **d. Record History**
     double current_Nm = 0.0, current_Nf = 0.0, current_Total_Matings = 0.0;
     for (int j = 0; j < num_block_types; ++j) {
       if (Current_Counts[j] > 0) {
         current_Nm += Males[j] * Current_Counts[j];
         current_Nf += Females[j] * Current_Counts[j];
         current_Total_Matings += Males[j] * Females[j] * Current_Counts[j];
       }
     }
     double current_Np = current_Nm + current_Nf;
     double current_sr = (current_Nf > 0) ? (current_Nm / current_Nf) : R_PosInf;
     double current_mean_mates = (current_Np > 0) ? (current_Total_Matings / current_Np) : 0.0;

     history_log_prob.push_back(Current_Log_Prob);
     history_Np.push_back(current_Np);
     history_sr.push_back(current_sr);
     history_mean_mates.push_back(current_mean_mates);

     // **e. Record Sample**
     if (t >= burn_in && ((t + 1 - burn_in) % thinning_interval == 0)) {
       results_list.push_back(clone(Current_Counts));
     }

     p.increment();
   }

   // --- Prepare Output ---
   double acceptance_rate = (double)accepted_count / n_iter;
   Rcout << "MCMC finished. Overall Acceptance Rate: " << round(acceptance_rate * 10000.0) / 100.0 << "%" << std::endl;

   // Create history DataFrame for return
   DataFrame history_df = DataFrame::create(
     _["iteration"] = seq(1, n_iter),
     _["log_prob"] = wrap(history_log_prob),
     _["Np"] = wrap(history_Np),
     _["sr"] = wrap(history_sr),
     _["mean_mates"] = wrap(history_mean_mates)
   );

   // Return results
   return List::create(
     _["samples"] = results_list,
     _["history"] = history_df,
     _["acceptance_rate"] = acceptance_rate
   );
 }
