#' @title Co-Mating Pair Density in a Bipartite Mating Network
#'
#' @description Quantifies micro-topological overlap as the proportion of all
#'   possible same-sex parent pairs that share at least one mating partner. This
#'   is the one-mode projection density of the bipartite mating graph, and it
#'   tracks how quickly local mating ties disperse across independent family
#'   lines as rewiring progresses.
#'
#' @details
#' The male projection is \code{mat \%*\% t(mat)} and the female projection is
#' \code{t(mat) \%*\% mat}, in both cases with the diagonal zeroed. A same-sex
#' pair is co-mating if its projection entry exceeds zero. Densities are reported
#' against the number of possible pairs of that sex, so the statistic is
#' comparable across population scales despite the quadratic growth of the dyadic
#' space.
#'
#' mateR2 matrices are always rows = males, columns = females. The projections
#' are computed directly from the matrix, so this function carries no igraph
#' dependency.
#'
#' @param mat A binary mating matrix (rows = males, columns = females). Non-zero
#'   entries are treated as edges.
#'
#' @return A single-row \code{data.frame} with:
#'   \item{male_density}{Proportion of male-male pairs sharing >= 1 female mate.}
#'   \item{female_density}{Proportion of female-female pairs sharing >= 1 male mate.}
#'   \item{mean_density}{Unweighted mean of the two, preventing bias toward the
#'     majority sex under skewed sex ratios.}
#'
#' @seealso \code{\link{count_network_components}}
#'
#' @examples
#' # Two males mated to the same single female: the male pair overlaps, and
#' # there is no female pair at all.
#' comating_pair_density(matrix(c(1, 1), nrow = 2, ncol = 1))
#'
#' @export
comating_pair_density <- function(mat) {
  mat <- as.matrix(mat)
  mat_bin <- matrix(as.numeric(mat != 0), nrow = nrow(mat), ncol = ncol(mat))

  n_males   <- nrow(mat_bin)
  n_females <- ncol(mat_bin)

  proj_males <- mat_bin %*% t(mat_bin)
  diag(proj_males) <- 0
  proj_females <- t(mat_bin) %*% mat_bin
  diag(proj_females) <- 0

  possible_male_pairs   <- if (n_males > 1)   n_males * (n_males - 1) / 2     else NA_real_
  possible_female_pairs <- if (n_females > 1) n_females * (n_females - 1) / 2 else NA_real_

  male_density   <- (sum(proj_males > 0) / 2) / possible_male_pairs
  female_density <- (sum(proj_females > 0) / 2) / possible_female_pairs

  data.frame(
    male_density   = male_density,
    female_density = female_density,
    mean_density   = mean(c(male_density, female_density), na.rm = TRUE)
  )
}
