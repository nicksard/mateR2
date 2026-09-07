make_ensemble <- function() {
  m <- generate_map_table(Np_target = 60, sr_target = 2, mean_mates_target = 1.5,
        max_males_per_female = 4, max_females_per_male = 4,
        n_iter = 5e4, burn_in = 5e3, thin = 50, seed = 1,
        verbose = FALSE, show_progress = FALSE)
  sample_posterior_ensemble(m, target_Np = 60, target_SR = 2, target_MM = 1.5,
        burn_in = 5e3, thin = 50, n_samples = 3, max_error_pct = 0.15)
}

test_that("a NULL sample size returns the complete cohort", {
  set.seed(1)
  e <- make_ensemble()
  full <- simulate_pedigree_ensemble(e, juvenile_sample_size = NULL,
            keep_matrices = TRUE, min_fecundity = 5, max_fecundity = 20)
  for (ped in full) {
    expect_equal(nrow(ped), sum(attr(ped, "K")))   # every offspring, none dropped
  }
})

test_that("a numeric sample size still sub-samples", {
  set.seed(2)
  e <- make_ensemble()
  s <- simulate_pedigree_ensemble(e, juvenile_sample_size = 50,
         min_fecundity = 5, max_fecundity = 20)
  expect_equal(nrow(s[[1]]), 50)
})

test_that("keep_matrices attaches M and K, and K is the COLONY mating matrix", {
  set.seed(3)
  e <- make_ensemble()
  ped <- simulate_pedigree_ensemble(e, juvenile_sample_size = NULL,
           keep_matrices = TRUE, min_fecundity = 5, max_fecundity = 20)[[1]]
  M <- attr(ped, "M"); K <- attr(ped, "K")
  expect_true(is.matrix(M) && is.matrix(K))
  expect_identical(dim(M), dim(K))
  expect_true(all((K > 0) == (M == 1)))          # support preserved
  out <- write_colony_sim(K, tempfile(fileext = ".Par"), mode = "sibship",
                          n_loci = 8, write_key = FALSE)
  expect_equal(out$n_offspring, sum(K))
})

test_that("matrices are withheld unless asked for", {
  set.seed(4)
  e <- make_ensemble()
  ped <- simulate_pedigree_ensemble(e, juvenile_sample_size = 50,
           min_fecundity = 5, max_fecundity = 20)[[1]]
  expect_null(attr(ped, "K"))
})

test_that("the package no longer calls the rewiring Curveball", {
  # Section 2.3 cites Gotelli & Entsminger and Miklos & Podani for a pairwise
  # checkerboard swap. Curveball (Strona et al. 2014) is a different algorithm,
  # and the mismatch is the first thing a reviewer opening the repo would see.
  rd <- list.files(system.file("man", package = "mateR2"), full.names = TRUE)
  hits <- unlist(lapply(rd, function(f) grep("Curveball", readLines(f), value = TRUE)))
  expect_true(all(grepl("not the Curveball", hits)))
})
