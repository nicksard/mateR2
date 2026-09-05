rt <- function(tbl) {
  m   <- mp_table_to_matrix(tbl)
  got <- mat2mp_sum_table(m)
  exp <- tbl[order(tbl$Males, tbl$Females), c("Males", "Females", "Count")]
  rownames(exp) <- NULL
  expect_equal(as.data.frame(got[, 1:3]), exp, ignore_attr = TRUE)
}

test_that("the block table round-trips at every sex ratio", {
  rt(data.frame(Males = c(1, 2), Females = c(1, 2), Count = c(2, 1)))  # square
  rt(data.frame(Males = 2,       Females = 1,       Count = 4))        # SR 2
  rt(data.frame(Males = 4,       Females = 1,       Count = 3))        # SR 4
  rt(data.frame(Males = c(2, 1), Females = c(1, 1), Count = c(3, 2)))  # SR 1.6
  rt(data.frame(Males = c(1, 2, 3), Females = c(1, 3, 2), Count = c(10, 4, 2)))
})

test_that("non-square matrices are handled rather than recycled or refused", {
  # Pairing rowSums with colSums positionally recycles silently whenever N_M is
  # a multiple of N_F, which is every SR = 2 and SR = 4 cell of the grid, and
  # errors otherwise. Blocks come from connectivity instead.
  for (d in list(c(8, 4), c(16, 4), c(6, 4), c(1067, 533))) {
    m <- matrix(0L, d[1], d[2])
    m[cbind(seq_len(d[1]), rep(seq_len(d[2]), length.out = d[1]))] <- 1L
    expect_s3_class(mat2mp_sum_table(m), "data.frame")
  }
})

test_that("blocks are counted, not individuals", {
  # Two 1:1 pairs plus one 2:2 block: the 2:2 block is one block containing two
  # males, and must not be counted twice.
  m <- mp_table_to_matrix(data.frame(Males = c(1, 2), Females = c(1, 2),
                                     Count = c(2, 1)))
  out <- mat2mp_sum_table(m)
  expect_equal(out$Count[out$Males == 2 & out$Females == 2], 1L)
  expect_equal(out$Count[out$Males == 1 & out$Females == 1], 2L)
})

test_that("degree signature alone cannot identify a block", {
  # Both matrices give every individual exactly one mate, but one is two 1:1
  # pairs and the other is a single 2:2 block with different connectivity.
  a <- mp_table_to_matrix(data.frame(Males = 1, Females = 1, Count = 2))
  b <- matrix(c(1L, 0L, 0L, 1L), 2, 2)[, 2:1]
  expect_equal(sort(rowSums(a)), sort(rowSums(b)))
  expect_equal(mat2mp_sum_table(a)$Count, mat2mp_sum_table(b)$Count)
})

test_that("unmated individuals are reported and excluded", {
  m <- mp_table_to_matrix(data.frame(Males = 1, Females = 1, Count = 3))
  m <- cbind(rbind(m, matrix(0L, 2, ncol(m))), matrix(0L, nrow(m) + 2, 1))
  out <- suppressMessages(mat2mp_sum_table(m))
  expect_equal(attr(out, "n_unmated_males"), 2L)
  expect_equal(attr(out, "n_unmated_females"), 1L)
  expect_equal(out$Count, 3L)
  kept <- suppressMessages(mat2mp_sum_table(m, drop_unmated = FALSE))
  expect_true(any(kept$Females == 0))
})

test_that("rewired matrices are flagged as no longer block-diagonal", {
  set.seed(1)
  b <- mp_table_to_matrix(data.frame(Males = c(2, 1), Females = c(2, 1),
                                     Count = c(6, 4)))
  expect_true(attr(mat2mp_sum_table(b), "complete_blocks"))
  mixed <- randomize_mating_structure(b, I = 0.5)
  expect_warning(out <- mat2mp_sum_table(mixed), "not complete bipartite")
  expect_false(attr(out, "complete_blocks"))
})

test_that("malformed input is refused", {
  expect_error(mat2mp_sum_table(matrix(0.5, 2, 2)), "0s and 1s")
  expect_error(mat2mp_sum_table(matrix(NA_integer_, 2, 2)), "0s and 1s")
  expect_error(mat2mp_sum_table(matrix(0L, 0, 2)), "zero dimension")
  expect_warning(mat2mp_sum_table(matrix(0L, 3, 3)), "no mating links")
})
