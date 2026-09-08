fx  <- function(f) testthat::test_path("fixtures", f)
sib <- function() read_colony_dat(fx("sim_sibship.DAT"))
par_ <- function() read_colony_dat(fx("sim_parentage.DAT"))

test_that("one row per locus, carrying the simulated rates when markers given", {
  d <- sib()
  r <- colony_locus_stats(d$offspring, d$markers)
  expect_equal(nrow(r$loci), nrow(d$markers))
  expect_equal(r$loci$marker, d$markers$marker)
  expect_equal(r$loci$dropout, d$markers$dropout)
  expect_equal(r$loci$error, d$markers$error)
})

test_that("markers are optional; names fall back to the columns", {
  d <- sib()
  r <- colony_locus_stats(d$offspring)
  expect_equal(r$loci$marker, d$markers$marker)   # stripped from Mk1-1 / Mk1-2
  expect_false("error" %in% names(r$loci))
})

test_that("a group label is added when asked for", {
  d <- sib()
  r <- colony_locus_stats(d$offspring, d$markers, group = "offspring")
  expect_equal(names(r$loci)[1], "group")
  expect_true(all(r$loci$group == "offspring"))
  expect_true(all(r$freqs$group == "offspring"))
})

test_that("allele frequencies sum to one within every locus", {
  r <- colony_locus_stats(sib()$offspring)
  s <- tapply(r$freqs$freq, r$freqs$marker, sum)
  expect_true(all(abs(s - 1) < 1e-9))
  n <- tapply(r$freqs$count, r$freqs$marker, sum)
  expect_equal(as.integer(n), 2L * r$loci$n_typed[match(names(n), r$loci$marker)])
})

test_that("only observed alleles appear, and n_alleles agrees", {
  r <- colony_locus_stats(sib()$offspring)
  expect_true(all(r$freqs$count > 0))
  obs <- tapply(r$freqs$allele, r$freqs$marker, length)
  expect_equal(as.integer(obs), r$loci$n_alleles[match(names(obs), r$loci$marker)])
})

test_that("missing genotypes are excluded, not counted as homozygotes", {
  d <- sib()
  n <- nrow(d$offspring)
  d$offspring[1:5, 2] <- 0L                 # zero one allele at locus 1
  r <- colony_locus_stats(d$offspring, d$markers)
  expect_equal(r$loci$n_missing[1], 5L)
  expect_equal(r$loci$n_typed[1], n - 5L)
  expect_equal(r$loci$missing_rate[1], 5 / n)
  expect_equal(r$loci$n_missing[2], 0L)     # other loci untouched
  expect_false(any(r$freqs$allele == 0))
})

test_that("heterozygosity is undefined rather than zero below two individuals", {
  d <- sib()
  d$offspring[-1, -1] <- 0L                 # one typed individual only
  r <- colony_locus_stats(d$offspring)
  expect_true(all(r$loci$n_typed == 1))
  expect_true(all(is.na(r$loci$He)))
})

test_that("Ho and He are bounded on real output", {
  r <- colony_locus_stats(par_()$offspring)
  expect_true(all(r$loci$Ho >= 0 & r$loci$Ho <= 1))
  expect_true(all(r$loci$He >= 0 & r$loci$He <= 1))
  expect_true(all(r$loci$n_alleles >= 1))
})

test_that("parents can be summarised alongside offspring", {
  d <- par_()
  all_ <- rbind(
    colony_locus_stats(d$offspring,         d$markers, group = "offspring")$loci,
    colony_locus_stats(d$candidate_males,   d$markers, group = "male")$loci,
    colony_locus_stats(d$candidate_females, d$markers, group = "female")$loci)
  expect_equal(nrow(all_), 3L * nrow(d$markers))
  expect_setequal(unique(all_$group), c("offspring", "male", "female"))
})

test_that("empty and malformed input are handled", {
  d <- sib()
  e <- colony_locus_stats(d$candidate_males, d$markers)   # sibship-only: 0 rows
  expect_equal(nrow(e$loci), 0)
  expect_error(colony_locus_stats(d), "not the whole file object")
  expect_error(colony_locus_stats(d$offspring[, 1:2]), "two columns per locus")
  expect_error(colony_locus_stats(d$offspring, d$markers[1:3, ]), "describes 3 loci")
})

test_that("freqs = FALSE returns the summary alone", {
  r <- colony_locus_stats(sib()$offspring, freqs = FALSE)
  expect_null(r$freqs)
  expect_equal(nrow(r$loci), 12)
})
