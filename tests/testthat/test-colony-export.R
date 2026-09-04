## Reduce a COLONY parameter file to the token stream the program actually
## reads: strip "!" notes, collapse whitespace, normalise numerics so that
## 0 and 0.0 compare equal.
tokenise <- function(lines) {
  x <- sub("!.*$", "", lines)
  x <- unlist(strsplit(trimws(x), "[ \t,]+"))
  x <- x[nzchar(x)]
  # Format element-wise: format() on a vector pads to a common number of
  # decimals, which would turn 0.05 into 0.0500 and break the comparison.
  vapply(x, function(tok) {
    n <- suppressWarnings(as.numeric(tok))
    if (is.na(n)) tok else format(n, trim = TRUE)
  }, character(1), USE.NAMES = FALSE)
}

demo_K <- function() {
  M <- matrix(0L, 8, 4,
              dimnames = list(sprintf("M%03d", 1:8), sprintf("F%03d", 1:4)))
  M[cbind(1:8, c(1, 1, 2, 2, 3, 3, 4, 4))] <- 1L
  M * rep(c(30, 20, 40, 25), each = 8)
}

test_that("the writer reproduces the COLONY user guide's worked example", {
  # Section 6.2 of the guide, transcribed. If the token stream matches, the
  # line order and content of the emitted file are correct by construction.
  guide <- c(
    "SimuTest !Project name", "2 !Number of replicates", "1 !method",
    "0 !precision", "2 0 !Dioecious, selfing rate", "1 !Number of matrices",
    "5 5 !#dads & mums",
    "1 0 0 0 0", "0 2 0 0 0", "0 0 4 0 0", "0 0 0 8 0", "0 0 0 0 16",
    "0 0 0 0 0", "0 0 0 0 0", "0 0 0 0 0", "0 0 0 0 0", "0 0 0 0 1",
    "0 0 0 0 0", "0 1 0 0 0", "0 0 0 0 0", "0 0 0 0 0", "0 0 0 0 0",
    "0 0 0 0 0", "0 0 0 0 0", "0 0 1 0 0", "0 0 0 0 0", "0 0 0 0 0",
    "1 0 0 0 0 !Dad included (1) or excluded (0)",
    "0 0 0 0 0 !mum included (1) or excluded (0)",
    "10 10 !# unrelated candidates", "0.50 0.50 !assumed prbs",
    "10 !Number of Loci", "0.05 !Prob. missing",
    ".01 .02 0 0 0 0 0 0 0 0 !Drop rate",
    ".01 0 .01 0 0 0 0 0 0 0 !OtherErrorRate",
    "0 0 0 0 0 0 0 0 0 0 !Marker types", "9 9 9 9 9 9 9 9 9 9 !# alleles",
    "0 !allele freq. distr.", "2 !ploidy", "1 1 !Monogamy/Polygamy",
    "1234 !Seed", "0 0.0 0.0 !sibship prior", "0 !Clone", "1 !Scale sibship",
    "0 !known allele frequency", "0 !updating allele freq", "1 !#replicate runs",
    "0 !run length", "-1 !Map length", "0.0 !Inbreeding coefficient",
    "0 !allowing inbreeding", "0 !DOS/GUI", "100000 !monitor")

  K  <- diag(c(1, 2, 4, 8, 16))
  kb <- matrix(0, 5, 5); kb[5, 5] <- 1
  kd <- matrix(0, 5, 5); kd[2, 2] <- 1
  km <- matrix(0, 5, 5); km[3, 3] <- 1

  out <- write_colony_sim(
    K, path = tempfile(fileext = ".Par"), mode = "parentage",
    project = "SimuTest", n_replicates = 2, likelihood = 1, precision = 0,
    known_both = kb, known_dad = kd, known_mum = km,
    sampled_dads = c(1, 0, 0, 0, 0), sampled_mums = rep(0, 5),
    n_candidate_males = 10, n_candidate_females = 10,
    prob_dad = 0.50, prob_mom = 0.50,
    n_loci = 10, missing_rate = 0.05,
    dropout_rate = c(.01, .02, rep(0, 8)),
    error_rate   = c(.01, 0, .01, rep(0, 7)),
    marker_type = 0, n_alleles = 9, allele_dist = "uniform",
    ploidy = 2, male_monogamous = TRUE, female_monogamous = TRUE,
    seed = 1234, sibship_prior = 0,
    paternal_sib_size = 0, maternal_sib_size = 0,
    clone = FALSE, scale_sibship = TRUE, n_runs = 1, run_length = 0,
    map_length = -1, inbreeding_coef = 0, inbreeding = FALSE,
    gui_mode = FALSE, monitor_interval = 100000, write_key = FALSE)

  expect_identical(tokenise(out$lines), tokenise(guide))
})

test_that("sibship mode emits no candidate parents", {
  K <- demo_K()
  out <- write_colony_sim(K, tempfile(fileext = ".Par"), mode = "sibship",
                          write_key = FALSE)
  # 4 matrices of nrow(K) rows, then the two parent-flag rows.
  hdr   <- grep("dads & mums", out$lines)
  flags <- out$lines[hdr + 2 + 4 * (nrow(K) + 1) + 0:1]
  expect_equal(tokenise(flags[1]), rep("0", nrow(K)))
  expect_equal(tokenise(flags[2]), rep("0", ncol(K)))
  expect_equal(tokenise(grep("unrelated candidate", out$lines, value = TRUE)),
               c("0", "0"))
})

test_that("parentage mode marks every true parent as sampled by default", {
  K <- demo_K()
  out <- write_colony_sim(K, tempfile(fileext = ".Par"), mode = "parentage",
                          n_candidate_males = 20, n_candidate_females = 20,
                          write_key = FALSE)
  hdr   <- grep("dads & mums", out$lines)
  flags <- out$lines[hdr + 2 + 4 * (nrow(K) + 1) + 0:1]
  expect_equal(tokenise(flags[1]), rep("1", nrow(K)))
  expect_equal(tokenise(flags[2]), rep("1", ncol(K)))
})

test_that("sibship mode refuses anything that would place a parent in the sample", {
  K <- demo_K()
  expect_error(write_colony_sim(K, tempfile(), mode = "sibship",
                                n_candidate_males = 5), "candidate parent")
  expect_error(write_colony_sim(K, tempfile(), mode = "sibship",
                                sampled_dads = rep(1, 8)), "sampled")
  expect_error(write_colony_sim(K, tempfile(), mode = "sibship",
                                known_both = (K > 0) * 1), "Known parentage")
})

test_that("the mating matrix round-trips through the file exactly", {
  K   <- demo_K()
  p   <- tempfile(fileext = ".Par")
  out <- write_colony_sim(K, p, mode = "sibship", write_key = FALSE)
  L   <- readLines(p)
  expect_identical(out$lines, L)             # nothing lost or added on write
  hdr <- grep("dads & mums", L)
  back <- as.matrix(utils::read.table(text = paste(L[hdr + 2:(nrow(K) + 1)],
                                                   collapse = "\n")))
  dimnames(back) <- dimnames(K)
  expect_equal(back, K)
})

test_that("a single write leaves no stray or doubled blank lines", {
  K <- demo_K()
  p <- tempfile(fileext = ".Par")
  write_colony_sim(K, p, mode = "sibship", write_key = FALSE)
  L <- readLines(p)
  expect_equal(sum(L[-1] == "" & utils::head(L, -1) == ""), 0)
  expect_true(all(nzchar(trimws(L[nzchar(L)]))))
})

test_that("the identifier key maps matrix position to mateR2 identifiers", {
  K   <- demo_K()
  p   <- tempfile(fileext = ".Par")
  out <- write_colony_sim(K, p, mode = "sibship", write_key = TRUE)
  key <- utils::read.csv(out$key_path, stringsAsFactors = FALSE)
  expect_equal(nrow(key), nrow(K) + ncol(K))
  expect_equal(key$mateR2_id[key$sex == "male"],   rownames(K))
  expect_equal(key$mateR2_id[key$sex == "female"], colnames(K))
  expect_equal(key$offspring[key$sex == "male"],   unname(rowSums(K)))
  expect_equal(key$offspring[key$sex == "female"], unname(colSums(K)))
})

test_that("invalid mating matrices are refused with an informative error", {
  K <- demo_K()
  expect_error(write_colony_sim(K + 0.5, tempfile()), "whole numbers")
  expect_error(write_colony_sim({z <- K; z[1, 1] <- -1; z}, tempfile()),
               "negative")
  expect_error(write_colony_sim(K * 0, tempfile()), "no offspring")
  expect_error(write_colony_sim({z <- K; z[1, 1] <- NA; z}, tempfile()),
               "missing values")
})

test_that("marker and project arguments are validated", {
  K <- demo_K()
  expect_error(write_colony_sim(K, tempfile(), project = "my sim.v2"),
               "letters and numbers")
  expect_error(write_colony_sim(K, tempfile(), n_loci = 10,
                                n_alleles = c(4, 4)), "one value per locus")
  expect_error(write_colony_sim(K, tempfile(), n_loci = 2, n_alleles = 5,
                                marker_type = 1), "Dominant markers")
  expect_error(write_colony_sim(K, tempfile(), mode = "parentage",
                                known_both = K * 0 + 1e6), "known-parentage")
})

test_that("empirical allele frequencies are validated and written per locus", {
  K <- demo_K()
  expect_error(write_colony_sim(K, tempfile(), allele_dist = "empirical"),
               "requires 'allele_freqs'")
  expect_error(write_colony_sim(K, tempfile(), n_loci = 2, n_alleles = 2,
                                allele_dist = "empirical",
                                allele_freqs = list(c(.5, .5), c(.2, .2))),
               "sum to")
  expect_error(write_colony_sim(K, tempfile(), n_loci = 2, n_alleles = c(2, 3),
                                allele_dist = "empirical",
                                allele_freqs = list(c(.5, .5), c(.5, .5))),
               "n_alleles says")

  out <- write_colony_sim(K, tempfile(fileext = ".Par"), n_loci = 2,
                          n_alleles = c(2, 3), allele_dist = "empirical",
                          allele_freqs = list(c(.1, .9), c(.05, .1405, .8095)),
                          write_key = FALSE)
  i <- grep("allele freq. distr", out$lines)
  expect_equal(tokenise(out$lines[i]), "3")
  expect_equal(tokenise(out$lines[i + 1]), c("0.1", "0.9"))
  expect_equal(tokenise(out$lines[i + 2]), c("0.05", "0.1405", "0.8095"))
})

test_that("wrapping splits rows without changing the token stream", {
  K <- demo_K()
  a <- write_colony_sim(K, tempfile(fileext = ".Par"), mode = "sibship",
                        seed = 7, write_key = FALSE)
  b <- write_colony_sim(K, tempfile(fileext = ".Par"), mode = "sibship",
                        seed = 7, write_key = FALSE, wrap_at = 2)
  expect_gt(length(b$lines), length(a$lines))
  expect_identical(tokenise(a$lines), tokenise(b$lines))
})

test_that("an existing file is not clobbered unless asked", {
  K <- demo_K()
  p <- tempfile(fileext = ".Par")
  write_colony_sim(K, p, mode = "sibship", write_key = FALSE)
  expect_error(write_colony_sim(K, p, mode = "sibship"), "already exists")
  expect_silent(suppressMessages(
    write_colony_sim(K, p, mode = "sibship", write_key = FALSE,
                     overwrite = TRUE)))
})

test_that("the seed obeys set.seed() when not supplied", {
  K <- demo_K()
  grab <- function() {
    out <- write_colony_sim(K, tempfile(fileext = ".Par"), mode = "sibship",
                            write_key = FALSE)
    tokenise(grep("Seed for random", out$lines, value = TRUE))
  }
  set.seed(99); a <- grab()
  set.seed(99); b <- grab()
  expect_identical(a, b)
})
