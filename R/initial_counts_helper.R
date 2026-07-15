#' @title Create Initial Counts for MCMC Sampler
#' @description This function generates an initial counts vector to serve as a
#'   "warm start" for the MCMC sampler. It dynamically selects the mating
#'   block configuration that best approximates the target sex ratio and
#'   populates the initial state to meet the target population size.
#'   This is a crucial first step for the Metropolis algorithm.
#' @param config_info The mating configuration table, typically created by
#'   `create_config_info()`.
#' @param Np_target The target total number of parents (males + females) for
#'   the simulation.
#' @param sr_target The target sex ratio (males / females). Defaults to 1.0.
#' @return A numeric vector representing the initial count of each mating
#'   configuration type.
#' @examples
#' # Create a sample config table and an initial counts vector
#' # config <- create_config_info(max_males_per_female = 3, max_females_per_male = 3)
#' # initial_counts <- create_initial_counts(config, Np_target = 100, sr_target = 1.0)
#' # print(initial_counts)
#' @export
create_initial_counts <- function(config_info, Np_target, sr_target = 1.0) {

  num_block_types <- nrow(config_info)
  initial_counts <- numeric(num_block_types)

  # Handle potential capitalization differences in package versions safely
  m_col <- intersect(c("Males", "males"), colnames(config_info))[1]
  f_col <- intersect(c("Females", "females"), colnames(config_info))[1]

  if (is.na(m_col) || is.na(f_col)) {
    stop("config_info must contain 'Males' and 'Females' columns.")
  }

  males <- config_info[[m_col]]
  females <- config_info[[f_col]]

  # 1. Calculate the sex ratio and total size of every block
  block_ratios <- males / females
  block_size <- males + females

  # 2. Create a searchable index dataframe
  calc_df <- data.frame(
    index = 1:num_block_types,
    ratio_diff = abs(block_ratios - sr_target),
    complexity = block_size
  )

  # 3. Primary sort: Closest to target SR
  #    Secondary sort (tie-breaker): Smallest block size (lowest complexity)
  calc_df <- calc_df[order(calc_df$ratio_diff, calc_df$complexity), ]

  # 4. Extract the winning block
  best_idx <- calc_df$index[1]
  best_males <- males[best_idx]
  best_females <- females[best_idx]
  best_block_size <- block_size[best_idx]
  best_ratio <- best_males / best_females

  # 5. Notify user if the warm start isn't a perfect match due to max_mates constraints
  if (abs(best_ratio - sr_target) > 0.05) {
    message(sprintf(
      "Note: The closest available block ratio (%d:%d = %.2f) differs from your target SR (%.2f). MCMC will adjust this organically.",
      best_males, best_females, best_ratio, sr_target
    ))
  }

  # 6. Calculate how many of these blocks are needed to reach Np_target
  count_best <- round(Np_target / best_block_size)
  if (count_best < 0) count_best <- 0

  initial_counts[best_idx] <- count_best

  return(initial_counts)
}
