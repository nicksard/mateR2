test_that("the MCMC sampler is reproducible from a seed", {
  args <- list(Np_target = 60, sr_target = 1, mean_mates_target = 1.5,
               max_males_per_female = 3, max_females_per_male = 3,
               np_weight = 50, sr_weight = 50, mm_weight = 50,
               n_iter = 5000, burn_in = 500, thin = 10, seed = 123)
  h1 <- do.call(generate_map_table, args)$mcmc_output$history
  h2 <- do.call(generate_map_table, args)$mcmc_output$history
  expect_identical(h1$Np, h2$Np)
  expect_identical(h1$log_prob, h2$log_prob)
})

test_that("an upstream set.seed() also fixes the chain", {
  args <- list(Np_target = 60, sr_target = 1, mean_mates_target = 1.5,
               max_males_per_female = 3, max_females_per_male = 3,
               n_iter = 3000, burn_in = 300, thin = 10)
  set.seed(7); h1 <- do.call(generate_map_table, args)$mcmc_output$history
  set.seed(7); h2 <- do.call(generate_map_table, args)$mcmc_output$history
  expect_identical(h1$Np, h2$Np)
})

test_that("different seeds give different chains", {
  args <- list(Np_target = 60, sr_target = 1, mean_mates_target = 1.5,
               max_males_per_female = 3, max_females_per_male = 3,
               n_iter = 3000, burn_in = 300, thin = 10)
  h1 <- do.call(generate_map_table, c(args, seed = 1))$mcmc_output$history
  h2 <- do.call(generate_map_table, c(args, seed = 2))$mcmc_output$history
  expect_false(identical(h1$Np, h2$Np))
})

test_that("block counts never go negative", {
  res <- generate_map_table(Np_target = 40, sr_target = 1, mean_mates_target = 1,
                            max_males_per_female = 3, max_females_per_male = 3,
                            n_iter = 5000, burn_in = 500, thin = 10, seed = 3)
  expect_true(all(res$map_table$MAP_Count >= 0))
})

test_that("the graph-theoretic floor MM >= (SR + 1) / 2 is enforced", {
  expect_error(check_target_viability(2, 1.0, 10, 10), "BOUNDARY")
  expect_error(check_target_viability(4, 2.0, 10, 10), "BOUNDARY")
  expect_true(check_target_viability(2, 1.5, 10, 10))
  expect_true(check_target_viability(4, 2.5, 10, 10))
  expect_true(check_target_viability(1, 1.0, 10, 10))
})

test_that("generate_map_table refuses infeasible targets before sampling", {
  expect_error(
    generate_map_table(Np_target = 100, sr_target = 4, mean_mates_target = 1,
                       max_males_per_female = 10, max_females_per_male = 10,
                       n_iter = 1000, burn_in = 100, thin = 10),
    "BOUNDARY"
  )
})

test_that("biological caps are applied to the correct sex", {
  cfg <- create_config_info(max_males_per_female = 3, max_females_per_male = 2)
  # A block's females each mate with all its males, so block males <= 3
  expect_equal(max(cfg$Males), 3)
  # A block's males each mate with all its females, so block females <= 2
  expect_equal(max(cfg$Females), 2)
  expect_equal(nrow(cfg), 6)
})

test_that("run_mcmc_chains disperses starts and returns diagnostics", {
  skip_if_not_installed("coda")
  r <- run_mcmc_chains(Np_target = 100, sr_target = 2, mean_mates_target = 2,
                       max_males_per_female = 10, max_females_per_male = 10,
                       n_chains = 3, n_iter = 8000, burn_in = 1000, thin = 10,
                       np_weight = 50, sr_weight = 50, mm_weight = 50, seed = 1)
  expect_equal(nrow(r$diagnostics), 3L)
  expect_setequal(r$diagnostics$parameter, c("Np", "sr", "mean_mates"))
  expect_length(r$chains, 3L)
  expect_true(all(is.finite(r$diagnostics$ess)))

  starts <- vapply(r$chains, function(x) sum(x$mcmc_output$history$Np[1]), numeric(1))
  expect_gt(length(unique(starts)), 1L)   # chains did not all start identically
})
