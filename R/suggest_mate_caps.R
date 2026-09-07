#' @title Suggest Mating Caps for a Set of Demographic Targets
#'
#' @description Converts demographic targets into the two mating caps
#'   \code{\link{generate_map_table}} and \code{\link{run_mcmc_chains}} require,
#'   set as a multiple of the mating success the targets imply rather than as a
#'   flat number, and reports the size of the state space that results.
#'
#' @details
#' In a fully connected block of \emph{m} males and \emph{f} females each male
#' has \emph{f} mates and each female has \emph{m}, so the caps are bounded below
#' by the mean mating success the targets imply:
#'
#' \preformatted{
#' mean mates of a male    =  2 * MM / (SR + 1)
#' mean mates of a female  =  2 * MM * SR / (SR + 1)
#' }
#'
#' \code{headroom} is the multiple of those requirements to allow. A flat cap is
#' not a constant amount of headroom: a cap of 10 against a requirement of one
#' mate is ten times the requirement, while the same cap against a requirement of
#' 5.3 is under twice, so a single number applied across a demographic grid poses
#' a different search problem in every cell.
#'
#' \strong{Why the default is 2.} Sweeping the headroom multiple across six
#' mating systems shows demographic accuracy improving monotonically as the caps
#' widen, so there is no accuracy argument for a tight cap. What degrades at
#' large multiples is mixing, and only because the caps define the state space:
#' the number of block types is their product, so widening both caps enlarges the
#' space quadratically while the iteration count stays fixed. A multiple of 1.5
#' to 2 captures nearly all of the available accuracy and sits at the peak of the
#' mixing curve. The \dQuote{Choosing Mate Caps} vignette shows the sweep.
#'
#' \strong{On the iteration warning.} Supplying \code{n_iter} triggers a rough
#' check of iterations per block type. The threshold is a rule of thumb drawn
#' from a limited sweep and is not a reliable predictor on its own: cells with
#' comparable ratios have converged and failed to converge. Treat it as a prompt
#' to inspect the diagnostics, never as a substitute for them.
#'
#' @param Np_target,sr_target,mean_mates_target The demographic targets, as
#'   passed to \code{\link{generate_map_table}}. \code{Np_target} is used only to
#'   report a suggested tolerance weight and does not affect the caps.
#' @param headroom Multiple of the implied mating success to allow. Default 2.
#' @param n_iter Optional planned iteration count. When supplied, the result
#'   carries a note on whether the run length looks adequate for the state space.
#'
#' @return A list with \code{max_males_per_female} and
#'   \code{max_females_per_male} ready to pass on, plus \code{block_types} (the
#'   size of the state space), \code{required_male_mates},
#'   \code{required_female_mates}, \code{suggested_weight} (the linear rule
#'   \code{50 * Np_target / 100}) and, when \code{n_iter} is given,
#'   \code{iter_per_block_type}.
#'
#' @seealso \code{\link{check_target_viability}} for the feasibility floor,
#'   \code{\link{generate_map_table}} for the sampler.
#'
#' @examples
#' suggest_mate_caps(Np_target = 400, sr_target = 2, mean_mates_target = 2)
#'
#' # A demanding cell, with the planned run length checked against it
#' suggest_mate_caps(400, sr_target = 1, mean_mates_target = 4, n_iter = 1e6)
#'
#' @export
suggest_mate_caps <- function(Np_target, sr_target, mean_mates_target,
                              headroom = 2, n_iter = NULL) {

  for (nm in c("Np_target", "sr_target", "mean_mates_target", "headroom")) {
    v <- get(nm)
    if (length(v) != 1L || !is.numeric(v) || is.na(v) || v <= 0) {
      stop(sprintf("\n[mateR2 INPUT ERROR]\n'%s' must be a single positive number.", nm),
           call. = FALSE)
    }
  }
  if (headroom < 1) {
    stop(sprintf(
      "\n[mateR2 INPUT ERROR]\nheadroom = %g would set the caps below the mating success the targets require, which cannot be satisfied.\nUse headroom of at least 1; 1.5 to 2 is the efficient range.",
      headroom), call. = FALSE)
  }

  req_male   <- 2 * mean_mates_target / (sr_target + 1)
  req_female <- 2 * mean_mates_target * sr_target / (sr_target + 1)

  if (req_male < 1) {
    stop(sprintf(
      "\n[mateR2 MATHEMATICAL BOUNDARY ERROR]\nTo achieve a global Mean Mates of %.2f at a Sex Ratio of %.2f, males must average %.2f mates.\nHowever, in a pedigree reconstruction, all successful males must have at least 1.0 mate.\n(Equivalently, Mean Mates must be at least (SR + 1)/2 = %.2f.)\nPlease raise your Mean Mates target or adjust your Sex Ratio skew.",
      mean_mates_target, sr_target, req_male, (sr_target + 1) / 2), call. = FALSE)
  }

  # A male's mates are females, so his cap is max_females_per_male.
  cap_f_per_m <- max(1L, as.integer(ceiling(headroom * req_male)))
  cap_m_per_f <- max(1L, as.integer(ceiling(headroom * req_female)))
  block_types <- cap_m_per_f * cap_f_per_m

  out <- list(
    max_males_per_female = cap_m_per_f,
    max_females_per_male = cap_f_per_m,
    block_types          = block_types,
    required_male_mates  = req_male,
    required_female_mates = req_female,
    headroom             = headroom,
    suggested_weight     = 50 * Np_target / 100
  )

  if (block_types == 1L) {
    message("These targets admit only 1:1 blocks, so the state space holds a single configuration. The sampler will find it exactly, and R-hat and effective sample size will be undefined because the within-chain variance is zero. That is not a failure.")
  }

  if (!is.null(n_iter)) {
    per <- n_iter / block_types
    out$iter_per_block_type <- per
    if (per < 2000) {
      warning(sprintf(
        "%s iterations across %d block types is about %.0f per type. Runs at that ratio have failed to converge in testing; consider raising n_iter, or lowering headroom if the caps are not biologically motivated. Check R-hat either way.",
        format(n_iter, big.mark = ",", scientific = FALSE), block_types, per),
        call. = FALSE)
    }
  }
  out
}
