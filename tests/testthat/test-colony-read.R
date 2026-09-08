fx <- function(f) testthat::test_path("fixtures", f)

test_that("the marker block is recovered", {
  d <- read_colony_dat(fx("sim_sibship.DAT"))
  expect_s3_class(d, "mateR2_colony_dat")
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

test_that("sibship-only files carry no candidate parents", {
  d <- read_colony_dat(fx("sim_sibship.DAT"))
  expect_equal(nrow(d$candidate_males), 0)
  expect_equal(nrow(d$candidate_females), 0)
  # the empty frames still have the genotype columns, so rbind stays valid
  expect_equal(ncol(d$candidate_males), 1 + 2 * 12)
})

test_that("candidate parents are read when present", {
  d <- read_colony_dat(fx("sim_parentage.DAT"))
  expect_equal(nrow(d$candidate_males), 16)     # 6 true + 10 unrelated
  expect_equal(nrow(d$candidate_females), 14)   # 4 true + 10 unrelated
  # COLONY marks true parents by case: M1 against m7
  expect_true(all(grepl("^M", d$candidate_males[[1]][1:6])))
  expect_true(all(grepl("^m", d$candidate_males[[1]][7:16])))
})

test_that("a file that is not a COLONY2 dat is refused", {
  p <- tempfile(); writeLines(c("just", "some", "text"), p)
  expect_error(read_colony_dat(p), "does not look like")
  expect_error(read_colony_dat(tempfile()), "does not exist")
})
