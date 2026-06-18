#' @title Assign Offspring Counts to a Breeding Matrix
#' @description This function populates a binary breeding matrix with offspring counts
#'   based on a specified fitness distribution type. It accounts for rounding
#'   to ensure the total number of offspring per female is preserved.
#' @param mat The binary breeding matrix (Individual-by-Individual Matrix).
#' @param min.fert The minimum total fecundity (offspring) for a female.
#' @param max.fert The maximum total fecundity (offspring) for a female.
#' @param type The type of offspring distribution among mates. Options are
#'   "uniform" (equal distribution) or "decline" (declining distribution with
#'   each successive mate).
#' @return A matrix of the same dimensions as the input, with offspring counts
#'   assigned to each mate pair.
#' @examples
#' # This function requires a binary matrix.
#' #
#' # # Example for a uniform distribution
#' # example_mat <- matrix(c(1,1,0,0, 1,1,0,0, 0,0,1,1, 0,0,1,1), nrow = 4, ncol = 4)
#' # mat_with_offspring <- brd.mat.fitness(example_mat, 1000, 2000, "uniform")
#' @export
brd.mat.fitness <- function(mat, min.fert, max.fert, type = "uniform") {
  # Check if input parameters are valid
  if (min.fert <= 0 || max.fert <= 0) {
    stop("min.fert and max.fert must be positive values.")
  }
  if (min.fert > max.fert) {
    stop("min.fert cannot be greater than max.fert.")
  }
  if (!(type %in% c("uniform", "decline"))) {
    stop("Invalid type argument. Must be 'uniform' or 'decline'.")
  }

  if (type == "uniform") {
    # Equal offspring distribution among mates
    for (i in 1:ncol(mat)) {
      my.mates <- which(mat[, i] > 0)
      my.off <- sample(x = min.fert:max.fert, size = 1)

      if (length(my.mates) > 1) {
        # Calculate equal proportions for each mate
        mate.props <- rep(1 / length(my.mates), length(my.mates))

        # Allocate offspring proportionally (with adjustment for rounding)
        my.off1 <- round(mate.props * my.off)
        diff <- my.off - sum(my.off1)
        if (diff != 0) {
          indices_to_adjust <- sample(1:length(my.off1), abs(diff))
          my.off1[indices_to_adjust] <- my.off1[indices_to_adjust] + sign(diff)
        }

        # Assign offspring counts to mates
        mat[my.mates, i] <- my.off1
      } else if (length(my.mates) == 1) {
        mat[my.mates, i] <- my.off
      }
    }
    return(mat)
  }

  if (type == "decline") {
    # Declining offspring distribution among mates
    for (i in 1:ncol(mat)) {
      my.mates <- which(mat[, i] > 0)
      my.off <- sample(x = min.fert:max.fert, size = 1)

      if (length(my.mates) > 1) {
        # Randomize the order of mates to avoid bias in declining proportions
        my.mates_shuffled <- sample(my.mates)

        # Calculate declining proportions for each mate
        mates_count <- length(my.mates_shuffled)
        mate.props <- (mates_count:1) / (mates_count * (mates_count + 1) / 2)

        # Allocate offspring proportionally (with adjustment for rounding)
        my.off1 <- round(mate.props * my.off)
        diff <- my.off - sum(my.off1)
        if (diff != 0) {
          indices_to_adjust <- sample(1:length(my.off1), abs(diff))
          my.off1[indices_to_adjust] <- my.off1[indices_to_adjust] + sign(diff)
        }

        # Assign offspring counts back to the original matrix order
        mat[my.mates_shuffled, i] <- my.off1
      } else if (length(my.mates) == 1) {
        mat[my.mates, i] <- my.off
      }
    }
    return(mat)
  }
}
