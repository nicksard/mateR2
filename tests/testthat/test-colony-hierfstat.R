fx <- function(f) testthat::test_path("fixtures", f)

## sim_sibship: 6 dads x 4 mums, every dad mated once (paternal nesting degenerate)
## sim_polygamous: 12 dads x 6 mums, dads mate twice and mums four times
sib  <- function() read_colony_dat(fx("sim_sibship.DAT"),  id_key = fx("sim_sibship_id_key.csv"))
poly <- function() read_colony_dat(fx("sim_polygamous.DAT"), id_key = fx("sim_polygamous_id_key.csv"))

test_that("the encoding round-trips through hierfstat's own decoder", {
  skip_if_not_installed("hierfstat")
  d  <- sib()
  hf <- colony_to_hierfstat(d, level = "full_sib")
  back <- hierfstat::getal.b(hf[, -1])
  orig <- as.matrix(d$offspring[, -1])
  a1 <- orig[, seq(1, ncol(orig), 2)]
  a2 <- orig[, seq(2, ncol(orig), 2)]
  expect_equal(unname(back[, , 1]), unname(pmin(a1, a2)))
  expect_equal(unname(back[, , 2]), unname(pmax(a1, a2)))
})

test_that("the digit width is fixed across the table, not per locus", {
  # getal.b picks its modulo from the global maximum, so a one-digit locus
  # mis-decodes as soon as any allele anywhere reaches 10. Two digits minimum.
  d  <- sib()
  hf <- colony_to_hierfstat(d, level = "full_sib")
  expect_equal(attr(hf, "ncode"), 2L)
  expect_true(all(hf[, -1] >= 100, na.rm = TRUE))
})

test_that("the frame has the shape hierfstat expects", {
  d  <- sib()
  hf <- colony_to_hierfstat(d, level = "full_sib")
  expect_equal(names(hf)[1], "pop")
  expect_equal(ncol(hf), 1 + nrow(d$markers))
  expect_equal(names(hf)[-1], d$markers$marker)
  expect_true(is.integer(hf$pop))
  expect_equal(attr(hf, "pop_labels")[1], "M001xF001")
})

test_that("grouping levels give the expected number of groups", {
  d <- poly()
  expect_equal(length(unique(colony_to_hierfstat(d, "full_sib")$pop)), 24)
  expect_equal(length(unique(colony_to_hierfstat(d, "paternal")$pop)), 12)
  expect_equal(length(unique(colony_to_hierfstat(d, "maternal")$pop)), 6)
  expect_equal(length(unique(colony_to_hierfstat(d, "cohort")$pop)), 1)
})

test_that("missing alleles become NA rather than a plausible genotype", {
  d <- sib()
  d$offspring[1:3, 2] <- 0L          # COLONY's missing code
  hf <- colony_to_hierfstat(d, "full_sib")
  expect_true(all(is.na(hf[1:3, 2])))
  expect_false(any(is.na(hf[1:3, 3])))
})

test_that("dominant markers are refused", {
  d <- sib()
  d$markers$type[3] <- 1L
  expect_error(colony_to_hierfstat(d, "full_sib"), "codominant")
})

test_that("min_size drops small families and keeps rows aligned", {
  d <- sib()
  full <- colony_to_hierfstat(d, "full_sib")
  cut  <- colony_to_hierfstat(d, "full_sib", min_size = 7)
  expect_lt(nrow(cut), nrow(full))
  expect_true(all(table(cut$pop) >= 7))
  expect_error(colony_to_hierfstat(d, "full_sib", min_size = 1000), "min_size")
})

test_that("degenerate nesting is caught rather than left to return NaN", {
  # Every sire in sim_sibship mated once, so full-sib families are not nested
  # within paternal families and varcomp.glob() would silently give NaN.
  expect_error(colony_hierfstat_levels(sib(), outer = "paternal"),
               "exactly one full-sibling family")
  expect_silent(colony_hierfstat_levels(sib(), outer = "maternal"))
})

test_that("levels are well formed when both sexes mate multiply", {
  d <- poly()
  for (s in c("paternal", "maternal")) {
    lev <- colony_hierfstat_levels(d, outer = s)
    expect_equal(nrow(lev), nrow(d$offspring))
    expect_equal(names(lev), c("parent", "full_sib"))
    # every full-sib family sits inside exactly one parental family
    per <- tapply(as.character(lev$parent), lev$full_sib,
                  function(v) length(unique(v)))
    expect_true(all(per == 1))
  }
})

test_that("statistics run and are ordered as the pedigree implies", {
  skip_if_not_installed("hierfstat")
  d <- poly()
  fst <- function(lv) hierfstat::basic.stats(colony_to_hierfstat(d, lv))$overall[["Fst"]]
  # Full sibs share more than half sibs, so among-family Fst is larger.
  expect_gt(fst("full_sib"), fst("paternal"))
  expect_gt(fst("full_sib"), fst("maternal"))
  # Within a family only a few parental gametes segregate, so heterozygote
  # excess is expected and Fis should be negative.
  expect_lt(hierfstat::basic.stats(colony_to_hierfstat(d, "full_sib"))$overall[["Fis"]], 0)
})
