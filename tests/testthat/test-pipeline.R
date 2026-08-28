test_that("mp_table_to_matrix reproduces the marginals implied by the count vector", {
  tab <- data.frame(Males = c(1, 2, 3), Females = c(1, 3, 2), Count = c(5, 2, 1))
  M <- mp_table_to_matrix(tab)
  expect_equal(nrow(M), sum(tab$Males * tab$Count))
  expect_equal(ncol(M), sum(tab$Females * tab$Count))
  expect_equal(sum(M), sum(tab$Males * tab$Females * tab$Count))
})

test_that("individual IDs are unique, sex-prefixed and externally valid", {
  M <- mp_table_to_matrix(data.frame(Males = 2, Females = 2, Count = 20))
  ids <- c(rownames(M), colnames(M))
  expect_length(unique(ids), length(ids))                 # no duplicates
  expect_length(intersect(rownames(M), colnames(M)), 0)   # none shared between sexes
  expect_true(all(grepl("^[A-Za-z0-9]+$", ids)))          # alphanumeric only
  expect_true(all(nchar(ids) <= 20))                      # COLONY's limit
})

test_that("edge swapping preserves every individual's degree", {
  M <- mp_table_to_matrix(data.frame(Males = 2, Females = 2, Count = 20))
  set.seed(1)
  Mm <- randomize_mating_structure(M, I = 1.0)
  # Names travel with the individual, so degrees must match name for name
  expect_equal(rowSums(M)[rownames(Mm)], rowSums(Mm))
  expect_equal(colSums(M)[colnames(Mm)], colSums(Mm))
  expect_equal(sum(M), sum(Mm))
  expect_true(all(Mm %in% c(0L, 1L)))
})

test_that("I = 0 leaves the mating structure intact", {
  M <- mp_table_to_matrix(data.frame(Males = 2, Females = 2, Count = 20))
  set.seed(2)
  M0 <- randomize_mating_structure(M, I = 0)
  # Only row/column ordering changes; restore it and the matrix is identical
  expect_equal(M0[rownames(M), colnames(M)], M)
})

test_that("mixing intensity is bounded", {
  M <- mp_table_to_matrix(data.frame(Males = 2, Females = 2, Count = 20))
  expect_error(randomize_mating_structure(M, I = -0.1), "bounded")
  expect_error(randomize_mating_structure(M, I = 1.5), "bounded")
})

test_that("IDs survive the whole pipeline without collision or re-indexing", {
  set.seed(11)
  M <- randomize_mating_structure(
    mp_table_to_matrix(data.frame(Males = 2, Females = 2, Count = 20)), I = 0.25)
  K <- brd.mat.fitness(M, min.fert = 100, max.fert = 200, type = "uniform")
  expect_identical(dimnames(K), dimnames(M))

  ped_full <- mat2ped(K)
  expect_true(all(ped_full$dad %in% rownames(K)))
  expect_true(all(ped_full$mom %in% colnames(K)))
  expect_length(intersect(ped_full$dad, ped_full$mom), 0)

  ped_samp <- convert2ped(mat.sub.sample(K, num_offspring = 200))
  expect_true(all(ped_samp$dad %in% rownames(K)))
  expect_length(intersect(ped_samp$dad, ped_samp$mom), 0)
})

test_that("mat2ped honours supplied dimnames rather than re-indexing", {
  K <- matrix(c(2, 0, 0, 3), 2, 2,
              dimnames = list(c("SIRE1", "SIRE2"), c("DAM1", "DAM2")))
  ped <- mat2ped(K)
  expect_setequal(unique(ped$dad), c("SIRE1", "SIRE2"))
  expect_setequal(unique(ped$mom), c("DAM1", "DAM2"))
  expect_equal(nrow(ped), sum(K))
})

test_that("fecundity allocation preserves the binary support", {
  set.seed(5)
  M <- randomize_mating_structure(
    mp_table_to_matrix(data.frame(Males = 2, Females = 2, Count = 15)), I = 0.5)
  K <- brd.mat.fitness(M, min.fert = 1000, max.fert = 5000, type = "uniform")
  expect_true(all((K > 0) == (M == 1)))
})

test_that("kbar_f = SR * kbar_m holds exactly", {
  set.seed(6)
  M <- randomize_mating_structure(
    mp_table_to_matrix(data.frame(Males = 2, Females = 1, Count = 30)), I = 0.25)
  K <- brd.mat.fitness(M, min.fert = 1000, max.fert = 5000, type = "uniform")
  NM <- nrow(K); NF <- ncol(K); SR <- NM / NF
  expect_equal(SR * (sum(K) / NM), sum(K) / NF)
})
