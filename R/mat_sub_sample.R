#' @title Sub-Sample Offspring from a Breeding Matrix (Null Model)
#' @description This function serves as the Null Model for field sampling. It
#'   simulates the process of capturing a fixed number of offspring from a
#'   breeding matrix under the assumption of perfect panmixia (no
#'   spatialtemporal clustering). The probability of sampling an offspring is
#'   strictly proportional to the initial number of offspring produced by its
#'   specific mate pair (a multinomial draw with replacement). This represents
#'   the "best-case scenario" for sampling offspring for pedigree analysis and
#'   serves as a baseline against which empirical or clustered sampling methods
#'   can be tested.
#' @param mat The breeding matrix with offspring counts.
#' @param num_offspring The number of offspring to sub-sample.
#' @return A data frame in long format containing the sampled mate pairs,
#'   their original offspring count, and the new sampled count.
#' @export
mat.sub.sample <- function(mat, num_offspring) {
  # Validates inputs
  if (any(mat < 0) || any(mat %% 1 != 0)) {
    stop("Input matrix 'mat' must contain non-negative integers.")
  }
  if (num_offspring <= 0) {
    stop("'num_offspring' must be a positive integer.")
  }
  if (num_offspring > sum(mat)) {
    warning("Sampling more offspring than available. The result will contain duplicate samples.")
  }

  # Convert matrix to long format with mate pair IDs (using base R)
  ped1 <- as.data.frame(as.table(mat))
  colnames(ped1) <- c("dads", "moms", "off")
  ped1 <- ped1[ped1$off != 0, ]
  ped1$mp <- paste(ped1$dads, ped1$moms, sep = "_")

  # Calculate sampling probabilities for each mate pair (vectorized)
  ped1$probs <- ped1$off / sum(ped1$off)

  # Randomly sample mate pairs based on probabilities
  sampled_mate_pairs <- sample(x = ped1$mp, size = num_offspring, replace = TRUE, prob = ped1$probs)

  # Count occurrences of each sampled mate pair
  mate_pair_counts <- as.data.frame(table(sampled_mate_pairs), stringsAsFactors = FALSE)
  names(mate_pair_counts) <- c("mp", "off1")

  # Merge with original data and fill missing counts with zero
  ped1 <- merge(x = ped1, y = mate_pair_counts, by = "mp", all.x = TRUE)
  ped1$off1[is.na(ped1$off1)] <- 0
  ped1$dads <- as.character(ped1$dads)
  ped1$moms <- as.character(ped1$moms)

  #order columns
  ped1 <- ped1[order(ped1$probs), ]

  return(ped1)
}
