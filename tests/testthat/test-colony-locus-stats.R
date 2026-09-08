fx <- function(f) testthat::test_path("fixtures", f)
sib  <- function() read_colony_dat(fx("sim_sibship.DAT"))
par_ <- function() read_colony_dat(fx("sim_parentage.DAT"))

test_that("one row per locus, carrying the simulated error rates", {
  d <- sib()
  r <- colony_locus_stats(d, scope = "offspring")
  expect_equal(nrow(r$loci), nrow(d$markers))
  expect_equal(r$loci$marker, d$markers$marker)
  expect_equal(r$loci$dropout, d$markers$dropout)
  expect_equal(r$loci$error, d$markers$error)
  expect_true(all(r$loci$group == "offspring"))
})

test_that("allele frequencies sum to one within every locus", {
  r <- colony_locus_stats(sib(), scope = "offspring")
  s <- tapply(r$freqs$freq, r$freqs$marker, sum)
  expect_true(all(abs(s - 1) < 1e-9))
  # counts are two per genotyped individual
  n <- tapply(r$freqs$count, r$freqs$marker, sum)
  expect_equal(as.integer(n), 2L * r$loci$n_typed[match(names(n), r$loci$marker)])
})

test_that("only observed alleles appear, and n_alleles agrees with the table", {
  r <- colony_locus_stats(sib(), scope = "offspring")
  expect_true(all(r$freqs$count > 0))
  obs <- tapply(r$freqs$allele, r$freqs$marker, length)
  expect_equal(as.integer(obs), r$loci$n_alleles[match(names(obs), r$loci$marker)])
})

test_that("missing genotypes are excluded, not counted as homozygotes", {
  d <- sib()
  n <- nrow(d$offspring)
  d$offspring[1:5, 2] <- 0L                 # zero one allele at locus 1
  r <- colony_locus_stats(d, scope = "offspring")
  expect_equal(r$loci$n_missing[1], 5L)
  expect_equal(r$loci$n_typed[1], n - 5L)
  expect_equal(r$loci$missing_rate[1], 5 / n)
  expect_equal(r$loci$n_missing[2], 0L)     # other loci untouched
  # a zero allele must never appear in the frequency table
  expect_false(any(r$freqs$allele == 0))
})

test_that("heterozygosity is undefined rather than zero below two individuals", {
  d <- sib()
  d$offspring[-1, seq(2, ncol(d$offspring))] <- 0L   # only one typed individual
  r <- colony_locus_stats(d, scope = "offspring")
  expect_true(all(r$loci$n_typed == 1))
  expect_true(all(is.na(r$loci$He)))
})

test_that("Ho and He are bounded and sane on real output", {
  r <- colony_locus_stats(par_(), scope = "offspring")
  expect_true(all(r$loci$Ho >= 0 & r$loci$Ho <= 1))
  expect_true(all(r$loci$He >= 0 & r$loci$He <= 1))
  expect_true(all(r$loci$n_alleles >= 1))
})

test_that("both scope returns offspring and each parent sex separately", {
  r <- colony_locus_stats(par_(), scope = "both")
  expect_setequal(unique(r$loci$group), c("offspring", "male", "female"))
  expect_equal(nrow(r$loci), 3L * nrow(par_()$markers))
  expect_setequal(unique(r$freqs$group), c("offspring", "male", "female"))
})

test_that("a sibship-only file has no parents to summarise", {
  expect_warning(r <- colony_locus_stats(sib(), scope = "parents"),
                 "no candidate parents")
  expect_equal(nrow(r$loci), 0)
})

test_that("freqs = FALSE returns the summary alone", {
  r <- colony_locus_stats(sib(), scope = "offspring", freqs = FALSE)
  expect_null(r$freqs)
  expect_equal(nrow(r$loci), 12)
})

test_that("input that did not come from read_colony_dat is refused", {
  expect_error(colony_locus_stats(list(offspring = data.frame())),
               "read_colony_dat")
})
