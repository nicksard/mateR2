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
create_config_info <- function(max_males_per_female, max_females_per_male,
                               complexity = c("sum", "cyclomatic", "quadratic", "asymmetry")) {
  complexity <- match.arg(complexity)
  if (max_males_per_female < 1 || max_females_per_male < 1) {
    stop("Maximum mates must be at least 1.")
  }
  config_info <- expand.grid(
    Males   = 1:max_males_per_female,
    Females = 1:max_females_per_male
  )
  config_info$Block <- paste(config_info$Males, config_info$Females, sep = ":")
  # "sum": delta_b = m + f - 2. Sums to N_P - 2 * n_blocks, so minimising it
  #   maximises the block count and favours many simple blocks.
  # "cyclomatic": delta_b = (m - 1)(f - 1), the number of independent cycles in
  #   the complete bipartite block. Zero for 1:1, 1:2 and 2:1; 80 for 11:9. This
  #   penalises communal (many-to-many) structure rather than block size.
  m <- config_info$Males; f <- config_info$Females
  # Note: N_P = sum(c_b (m + f)) and E = sum(c_b m f) are both pinned by the
  # demographic targets, so ANY delta linear or bilinear in m and f reduces to
  # +/- the block count and cannot express density. "sum" = N_P - 2 n_blocks;
  # "cyclomatic" = E - N_P + n_blocks. Only superlinear forms discriminate.
  config_info$Complexity_Diff <- switch(complexity,
    sum        = (m + f) - 2,
    cyclomatic = (m - 1) * (f - 1),
    quadratic  = ((m + f) - 2)^2,
    asymmetry  = abs(m - f)
  )
  config_info <- config_info[, c("Block", "Males", "Females", "Complexity_Diff")]
  config_info <- config_info[order(config_info$Males, config_info$Females), ]
  rownames(config_info) <- NULL
  return(config_info)
}
