#' @title Run Independent MCMC Chains and Assemble Convergence Diagnostics
#'
#' @description Runs several independent Metropolis-Hastings chains against the
#'   same demographic targets from dispersed starting states, and returns both the
#'   per-chain results and \code{coda} objects ready for Gelman-Rubin and
#'   effective sample size diagnostics. Running the chains through one function
#'   keeps seeding, dispersal, burn-in and thinning consistent across them.
#'
#' @details
#' Chain `i` uses seed `seed + i - 1`, so the whole ensemble is reproducible from
#' a single integer.
#'
#' \code{\link{create_initial_counts}} is deterministic, so every chain would
#' otherwise begin from the identical state and the Gelman-Rubin statistic would
#' be uninformative. Each chain therefore jitters the warm start: active block
#' counts are multiplied by a uniform draw between 0.6 and 1.4 and rounded, with
#' a floor of one. That disperses the starting states while keeping them near the
#' target manifold, so no chain begins somewhere infeasible.
#'
#' Traces are the post-burn-in portion of the recorded history, thinned by `thin`
#' and wrapped as \code{coda::mcmc} objects. Requires the coda package, which is
#' listed under Suggests.
#'
#' @param Np_target The target total number of parents (males + females).
#' @param sr_target The target sex ratio (males / females).
#' @param mean_mates_target The target mean number of mates per individual.
#' @param max_males_per_female The maximum number of males a female can mate with.
#' @param max_females_per_male The maximum number of females a male can mate with.
#' @param n_chains Integer. Number of independent chains (default 4).
#' @param seed Integer. Base seed; chain `i` uses `seed + i - 1` (default 1).
#' @param n_iter,burn_in,thin MCMC length, burn-in and thinning interval, passed
#'   to \code{\link{generate_map_table}}.
#' @param ... Further arguments passed to \code{\link{generate_map_table}}, such
#'   as `np_weight`, `sr_weight`, `mm_weight` and `decay_constant`.
#'
#' @return A list with:
#'   \item{chains}{List of the per-chain \code{\link{generate_map_table}} results.}
#'   \item{traces}{Named list of \code{coda::mcmc.list} objects, one per target
#'     parameter (`Np`, `sr`, `mean_mates`).}
#'   \item{diagnostics}{A \code{data.frame} of realised posterior means, percent
#'     error against target, Gelman-Rubin `rhat` and summed `ess`.}
#'
#' @examples
#' \dontrun{
#' res <- run_mcmc_chains(
#'   Np_target = 200, sr_target = 2, mean_mates_target = 2,
#'   max_males_per_female = 10, max_females_per_male = 10,
#'   n_chains = 4, n_iter = 50000, burn_in = 5000, thin = 10,
#'   np_weight = 70, sr_weight = 70, mm_weight = 70, decay_constant = -0.05
#' )
#' res$diagnostics
#' }
#'
#' @importFrom stats runif
#' @param trace_thin Thinning interval applied only to the returned diagnostic
#'   traces. Independent of \code{thin}, which governs the posterior sample
#'   pool; changing it leaves the chain itself bit-identical. Defaults to
#'   \code{thin}.
#' @param complexity Structural complexity metric assigned to each block type,
#'   one of "sum" (the default, m + f - 2), "cyclomatic" ((m - 1)(f - 1)),
#'   "quadratic" ((m + f - 2)^2) or "asymmetry" (|m - f|). Only the
#'   superlinear forms distinguish dense blocks from sparse ones of the same
#'   size.
#' @export
run_mcmc_chains <- function(Np_target, sr_target, mean_mates_target,
                            max_males_per_female, max_females_per_male,
                            n_chains = 4, seed = 1,
                            n_iter = 200000, burn_in = 20000, thin = 10,
                            trace_thin = NULL,
                            complexity = c("sum", "cyclomatic", "quadratic", "asymmetry"),
                            ...) {
  complexity <- match.arg(complexity)

  # `thin` governs the retained sample pool (MAP search + posterior ensemble).
  # `trace_thin` governs only the diagnostic traces. Unrelated choices; coupling
  # them makes ESS move when the chain itself has not changed.
  if (is.null(trace_thin)) trace_thin <- thin

  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for run_mcmc_chains(). ",
         "Install it with install.packages('coda').", call. = FALSE)
  }
  if (n_chains < 1) stop("'n_chains' must be at least 1.")

  # Validate once up front so an infeasible target fails immediately rather than
  # after the first chain has run.
  check_target_viability(sr_target, mean_mates_target,
                         max_males_per_female, max_females_per_male)

  config_info <- create_config_info(max_males_per_female, max_females_per_male,
                                    complexity = complexity)
  warm_start  <- create_initial_counts(config_info, Np_target, sr_target)

  chain_results <- vector("list", n_chains)

  for (i in seq_len(n_chains)) {
    chain_seed <- seed + i - 1

    # Disperse the starting state. create_initial_counts() is deterministic, so
    # without this every chain would start from the same point.
    set.seed(chain_seed)
    init <- warm_start
    active <- which(init > 0)
    init[active] <- pmax(1, round(init[active] *
                                    stats::runif(length(active), 0.6, 1.4)))

    chain_results[[i]] <- generate_map_table(
      Np_target = Np_target,
      sr_target = sr_target,
      mean_mates_target = mean_mates_target,
      max_males_per_female = max_males_per_female,
      max_females_per_male = max_females_per_male,
      n_iter = n_iter, burn_in = burn_in, thin = thin,
      complexity = complexity,
      initial_method = init,
      seed = chain_seed,
      ...
    )
  }
  names(chain_results) <- paste0("chain", seq_len(n_chains))

  # --- coda traces from the post-burn-in, thinned history ---
  params <- c(Np = "Np", sr = "sr", mean_mates = "mean_mates")
  traces <- lapply(params, function(p) {
    coda::mcmc.list(lapply(chain_results, function(res) {
      h    <- res$mcmc_output$history
      post <- h[h$iteration > burn_in, ]
      keep <- seq(1, nrow(post), by = trace_thin)
      coda::mcmc(post[[p]][keep], thin = trace_thin)
    }))
  })

  # --- Per-parameter diagnostics ---
  targets <- c(Np = Np_target, sr = sr_target, mean_mates = mean_mates_target)
  diagnostics <- do.call(rbind, lapply(names(params), function(p) {
    ml  <- traces[[p]]
    obs <- mean(vapply(ml, mean, numeric(1)))
    rhat <- if (n_chains > 1) {
      tryCatch(coda::gelman.diag(ml, autoburnin = FALSE,
                                 multivariate = FALSE)$psrf[1, 1],
               error = function(e) NA_real_)
    } else NA_real_
    ess <- tryCatch(sum(coda::effectiveSize(ml)), error = function(e) NA_real_)
    data.frame(
      parameter = p,
      target    = unname(targets[[p]]),
      observed  = obs,
      error_pct = abs(obs - targets[[p]]) / targets[[p]] * 100,
      rhat      = unname(rhat),
      ess       = unname(ess),
      n_chains  = n_chains,
      stringsAsFactors = FALSE
    )
  }))
  rownames(diagnostics) <- NULL

  list(chains = chain_results, traces = traces, diagnostics = diagnostics)
}
