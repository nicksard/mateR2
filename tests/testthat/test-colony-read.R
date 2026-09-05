## Fixtures are genuine Simu2.exe output, generated from a known mateR2 K.
fx <- function(f) testthat::test_path("fixtures", f)

sib_K <- function() {
  K <- matrix(0L, 6, 4,
              dimnames = list(sprintf("M%03d", 1:6), sprintf("F%03d", 1:4)))
  K[cbind(1:6, c(1, 1, 2, 2, 3, 4))] <- c(8L, 8L, 6L, 6L, 5L, 7L)
  K
}

test_that("the header and marker block are recovered", {
  d <- read_colony_dat(fx("sim_sibship.DAT"))
  expect_s3_class(d, "mateR2_colony_dat")
  expect_equal(d$params$n_offspring, 40)
  expect_equal(d$params$n_loci, 12)
  expect_equal(d$params$seed, 1234)
  expect_equal(nrow(d$markers), 12)
  expect_true(all(d$markers$type == 0))
  expect_equal(unique(d$markers$error), 0.001)
})

test_that("genotypes come back wide with two columns per locus", {
  d <- read_colony_dat(fx("sim_sibship.DAT"))
  expect_equal(nrow(d$offspring), 40)
  expect_equal(ncol(d$offspring), 1 + 2 * 12)
  expect_equal(names(d$offspring)[1:3], c("Offspring", "Mk1-1", "Mk1-2"))
  expect_true(all(vapply(d$offspring[-1], is.integer, logical(1))))
})

test_that("the true pedigree is recovered and matches the exported K", {
  d <- read_colony_dat(fx("sim_sibship.DAT"))
  expect_equal(nrow(d$truth), 40)
  expect_identical(colony_truth_matrix(d, 6, 4), unname(sib_K()))
})

test_that("the file's own configuration block agrees with the identifiers", {
  # Two independent routes to the same truth: COLONY's True Configuration
  # block, and the M<i>F<j>C<k> identifiers it generates.
  d <- read_colony_dat(fx("sim_sibship.DAT"))
  expect_equal(d$config$paternal_family, d$truth$dad_index)
  expect_equal(d$config$maternal_family, d$truth$mum_index)
})

test_that("the identifier key maps the pedigree back to mateR2 names", {
  d <- read_colony_dat(fx("sim_sibship.DAT"), id_key = fx("sim_sibship_id_key.csv"))
  expect_equal(d$truth$dad, sprintf("M%03d", d$truth$dad_index))
  expect_equal(d$truth$mum, sprintf("F%03d", d$truth$mum_index))
})

test_that("sibship-only files carry no candidate parents", {
  d <- read_colony_dat(fx("sim_sibship.DAT"))
  expect_equal(nrow(d$candidate_males), 0)
  expect_equal(nrow(d$candidate_females), 0)
})

test_that("candidate parents are read and true parents flagged by case", {
  d <- read_colony_dat(fx("sim_parentage.DAT"),
                       id_key = fx("sim_parentage_id_key.csv"))
  expect_equal(nrow(d$candidate_males), 16)     # 6 true + 10 unrelated
  expect_equal(nrow(d$candidate_females), 14)   # 4 true + 10 unrelated
  expect_equal(sum(d$candidate_males$true_parent), 6)
  expect_equal(sum(d$candidate_females$true_parent), 4)
  expect_equal(d$candidate_males$mateR2_id[1:6], sprintf("M%03d", 1:6))
  expect_true(all(is.na(d$candidate_males$mateR2_id[7:16])))
})

test_that("dat inputs are reshaped for an inference writer", {
  d <- read_colony_dat(fx("sim_parentage.DAT"),
                       id_key = fx("sim_parentage_id_key.csv"))
  x <- colony_dat_inputs(d, parents = "all")
  expect_equal(dim(x$markers), c(3, 10))        # type, dropout, error
  expect_equal(names(x$markers), d$markers$marker)
  expect_equal(nrow(x$dads), 16)
  expect_false("true_parent" %in% names(x$dads))
  expect_equal(x$dads[[1]][1], "M001")          # relabelled where known

  expect_equal(nrow(colony_dat_inputs(d, parents = "true")$dads), 6)
  expect_equal(nrow(colony_dat_inputs(d, parents = "none")$dads), 0)
})

test_that("sub-sampling keeps offspring and truth aligned", {
  d <- read_colony_dat(fx("sim_sibship.DAT"))
  set.seed(2)
  x <- colony_dat_inputs(d, n_offspring = 15, parents = "none")
  expect_equal(nrow(x$kids), 15)
  expect_equal(nrow(x$truth), 15)
  expect_equal(x$kids[[1]], x$truth$off)
  expect_error(colony_dat_inputs(d, n_offspring = 999), "file holds")
})

test_that("a file that is not a COLONY2 dat is refused", {
  p <- tempfile(); writeLines(c("just", "some", "text"), p)
  expect_error(read_colony_dat(p), "does not look like")
  expect_error(read_colony_dat(tempfile()), "does not exist")
})
