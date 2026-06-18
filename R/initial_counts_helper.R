#' @title Create Initial Counts for MCMC Sampler
#' @description This function generates an initial counts vector, primarily
#'   composed of 1:1 mating pairs, to serve as a starting point for the MCMC
#'   sampler. This is a crucial first step for the Metropolis algorithm.
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
  if (abs(sr_target - 1.0) > 1e-6) {
    warning("Creating initial counts with only 1:1 blocks assumes sr_target=1.0.")
  }
  num_block_types <- nrow(config_info)
  initial_counts <- numeric(num_block_types)
  index_11 <- which(config_info$Males == 1 & config_info$Females == 1)
  if (length(index_11) != 1) { stop("Could not find unique 1:1 block type.") }
  count_11 <- round(Np_target / 2)
  if (count_11 < 0) count_11 <- 0
  initial_counts[index_11] <- count_11
  start_Np <- sum(config_info$Males * initial_counts) + sum(config_info$Females * initial_counts)
  return(initial_counts)
}
