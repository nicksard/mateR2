#' @title Create a Mating Configuration Information Table
#' @description This function generates a data frame that serves as a lookup
#'   table for all possible breeding configurations (mate pairs) within a
#'   defined system. It is a necessary input for the MCMC sampler.
#' @param max_males_per_female The maximum number of males a single female can
#'   mate with.
#' @param max_females_per_male The maximum number of females a single male can
#'   mate with.
#' @return A data frame with columns: `Block` (e.g., "1:1"), `Males`, `Females`,
#'   and `Complexity_Diff`, which is used to penalize more complex mating
#'   configurations in the MCMC algorithm.
#' @examples
#' # Create a config table for a system where a single individual can have up to 3 mates
#' config_table <- create_config_info(max_males_per_female = 3, max_females_per_male = 3)
#' print(config_table)
#' @export
create_config_info <- function(max_males_per_female, max_females_per_male) {
  if (max_males_per_female < 1 || max_females_per_male < 1) {
    stop("Maximum mates must be at least 1.")
  }
  config_info <- expand.grid(
    Males = 1:max_females_per_male,
    Females = 1:max_males_per_female
  )
  config_info$Block <- paste(config_info$Males, config_info$Females, sep = ":")
  config_info$Complexity_Diff <- (config_info$Males + config_info$Females) - 2
  config_info <- config_info[, c("Block", "Males", "Females", "Complexity_Diff")]
  config_info <- config_info[order(config_info$Males, config_info$Females), ]
  rownames(config_info) <- NULL
  return(config_info)
}
