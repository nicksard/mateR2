test_that("sub-sampling is without replacement", {
  K <- matrix(c(30, 0, 10, 0, 40, 20), nrow = 2, byrow = TRUE)
  set.seed(1)
  for (n in c(5, 50, 100)) {
    s <- mat.sub.sample(K, num_offspring = n)
    expect_equal(sum(s$off1), n)
    expect_true(all(s$off1 <= s$off))   # never over-sampled
    expect_true(all(s$off1 >= 0))
  }
})

test_that("sub-sampling refuses to exceed the cohort", {
  K <- matrix(c(3, 0, 0, 2), 2, 2)
  expect_error(mat.sub.sample(K, num_offspring = 6), "SAMPLING ERROR")
  s <- mat.sub.sample(K, num_offspring = 5)
  expect_equal(sum(s$off1), 5)
})

test_that("sub-sampling retains pairs with no sampled offspring", {
  K <- matrix(c(1000, 0, 0, 1000), 2, 2)
  set.seed(4)
  s <- mat.sub.sample(K, num_offspring = 10)
  expect_equal(nrow(s), 2)   # both active pairs present even if one yields zero
})

test_that("sub-sampling gives every offspring the same inclusion probability", {
  # Two pairs at a 3:1 fecundity ratio; expected sampled counts follow that ratio
  K <- matrix(c(3000, 0, 0, 1000), 2, 2,
              dimnames = list(c("M1", "M2"), c("F1", "F2")))
  set.seed(9)
  draws <- replicate(200, {
    s <- mat.sub.sample(K, num_offspring = 400)
    s$off1[s$mp == "M1_F1"]
  })
  expect_equal(mean(draws), 300, tolerance = 0.02)
})

test_that("sib.stats per-offspring means agree with the dyad counts", {
  mat <- matrix(c(2, 0, 1, 0, 4, 3), nrow = 2, byrow = TRUE)
  s <- sib.stats(mat)
  n <- sum(mat)
  expect_equal(s$mean_FS,  2 * s$FS_pairs  / n)
  expect_equal(s$mean_MHS, 2 * s$MHS_pairs / n)
  expect_equal(s$mean_PHS, 2 * s$PHS_pairs / n)
})

test_that("sib.stats reports full-siblings per offspring", {
  # Four offspring in one cell: each has three full sibs and no half sibs
  s <- sib.stats(matrix(c(4, 0, 0, 0), 2, 2))
  expect_equal(s$mean_FS, 3)
  expect_equal(s$mean_MHS, 0)
  expect_equal(s$mean_PHS, 0)
})

test_that("sib.stats is invariant to labels and to row/column order", {
  set.seed(12)
  K <- brd.mat.fitness(
    randomize_mating_structure(
      mp_table_to_matrix(data.frame(Males = 2, Females = 1, Count = 20)), I = 0.25),
    min.fert = 100, max.fert = 500, type = "uniform")
  Ku <- K; dimnames(Ku) <- NULL
  expect_equal(sib.stats(K), sib.stats(Ku))
  Kp <- K[sample(nrow(K)), sample(ncol(K))]
  expect_equal(sib.stats(K), sib.stats(Kp))
})

test_that("parent.class.stats counts yield tiers per sex", {
  mat <- matrix(c(1, 0, 2, 0, 2, 0), nrow = 2, byrow = TRUE)
  p <- parent.class.stats(mat)
  expect_equal(p["maternal", "detected"], 3)
  expect_equal(p["maternal", "singletons"], 1)
  expect_equal(p["maternal", "doubletons"], 2)
  expect_equal(p["paternal", "detected"], 2)
})

test_that("mat.stats reports V_k and does not confuse it with k-bar", {
  K <- matrix(c(10, 0, 0, 30), 2, 2)
  s <- mat.stats(K)
  expect_equal(s$mean.male.rs, 20)
  expect_equal(s$var.male.rs, stats::var(c(10, 30)))
  expect_equal(s$offspring.per.parent, 10)   # 40 offspring / 4 parents
})

test_that("co-mating density is orientation-correct", {
  # Two males both mated to the same single female: the male pair overlaps, and
  # there is no female pair at all.
  d <- comating_pair_density(matrix(c(1, 1), nrow = 2, ncol = 1))
  expect_equal(d$male_density, 1)
  expect_true(is.na(d$female_density))
  expect_equal(d$mean_density, 1)
})

test_that("rewiring dissolves block structure", {
  skip_if_not_installed("igraph")
  M <- mp_table_to_matrix(data.frame(Males = 2, Females = 2, Count = 20))
  set.seed(3)
  expect_gt(count_network_components(M),
            count_network_components(randomize_mating_structure(M, I = 0.25)))
})
