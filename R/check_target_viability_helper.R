#' @title Check Demographic Target Viability
#' @description Evaluates whether a requested combination of Sex Ratio and Mean Mates
#'   is mathematically possible given the biological mating caps of the species.
#' @param sr_target The target operational sex ratio (Males / Females).
#' @param mean_mates_target The target global mean mates.
#' @param max_males_per_female Biological cap on how many males a female can mate with.
#' @param max_females_per_male Biological cap on how many females a male can mate with.
#' @return TRUE if viable. Stops execution with an informative error if impossible.
#' @export
check_target_viability <- function(sr_target, mean_mates_target,
                                   max_males_per_female, max_females_per_male) {

  # FIXED: Using the true Average of Averages Bipartite Algebra
  req_female_mates <- mean_mates_target * (2 * sr_target) / (sr_target + 1)
  req_male_mates   <- mean_mates_target * 2 / (sr_target + 1)

  # 1. LOWER BOUND CHECKS (Mathematical Parentage Floor)
  # Every individual in a parentage-reconstructed pedigree must have >= 1.0 mate.

  if (req_female_mates < 1.0) {
    stop(sprintf(
      "\n[mateR2 MATHEMATICAL BOUNDARY ERROR]\nTo achieve a global Mean Mates of %.2f at a Sex Ratio of %.2f, females must average %.2f mates.\nHowever, in a pedigree reconstruction, all successful females must have at least 1.0 mate.\nPlease raise your Mean Mates target or reduce your Sex Ratio skew.",
      mean_mates_target, sr_target, req_female_mates
    ), call. = FALSE)
  }

  if (req_male_mates < 1.0) {
    stop(sprintf(
      "\n[mateR2 MATHEMATICAL BOUNDARY ERROR]\nTo achieve a global Mean Mates of %.2f at a Sex Ratio of %.2f, males must average %.2f mates.\nHowever, in a pedigree reconstruction, all successful males must have at least 1.0 mate.\nPlease raise your Mean Mates target or adjust your Sex Ratio skew.",
      mean_mates_target, sr_target, req_male_mates
    ), call. = FALSE)
  }

  # 2. UPPER BOUND CHECKS (Biological Capacity Ceiling)

  # Check if females violate the biological cap
  if (req_female_mates > max_males_per_female) {
    stop(sprintf(
      "\n[mateR2 BIOLOGICAL PARADOX ERROR]\nTo achieve a global Mean Mates of %.2f at a Sex Ratio of %.2f, females must average %.2f mates.\nYour biological cap restricts females to a maximum of %d mates.\nPlease lower your Mean Mates target, reduce your Sex Ratio skew, or raise the maximum female capacity.",
      mean_mates_target, sr_target, req_female_mates, max_males_per_female
    ), call. = FALSE)
  }

  # Check if males violate the biological cap
  if (req_male_mates > max_females_per_male) {
    stop(sprintf(
      "\n[mateR2 BIOLOGICAL PARADOX ERROR]\nTo achieve a global Mean Mates of %.2f at a Sex Ratio of %.2f, males must average %.2f mates.\nYour biological cap restricts males to a maximum of %d mates.\nPlease lower your Mean Mates target, adjust your Sex Ratio skew, or raise the maximum male capacity.",
      mean_mates_target, sr_target, req_male_mates, max_females_per_male
    ), call. = FALSE)
  }

  return(TRUE)
}
