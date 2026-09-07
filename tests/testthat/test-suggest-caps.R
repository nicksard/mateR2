test_that("caps are never below the mating success the targets imply", {
  for (SR in c(1, 2, 4)) for (MM in c(1, 2, 4)) {
    if (MM < (SR + 1) / 2) next
    s <- suggest_mate_caps(400, SR, MM, headroom = 1)
    expect_gte(s$max_females_per_male, s$required_male_mates)
    expect_gte(s$max_males_per_female, s$required_female_mates)
  }
})

test_that("the caps map to the right sex", {
  # A male's mates are females, so his cap is max_females_per_male. At SR = 4
  # females mate far more often than males, so their cap must be the larger.
  s <- suggest_mate_caps(400, sr_target = 4, mean_mates_target = 4)
  expect_gt(s$required_female_mates, s$required_male_mates)
  expect_gt(s$max_males_per_female, s$max_females_per_male)
})

test_that("headroom widens the caps and block_types is their product", {
  a <- suggest_mate_caps(400, 2, 2, headroom = 1)
  b <- suggest_mate_caps(400, 2, 2, headroom = 4)
  expect_gt(b$max_males_per_female, a$max_males_per_female)
  expect_gt(b$max_females_per_male, a$max_females_per_male)
  for (s in list(a, b)) {
    expect_equal(s$block_types, s$max_males_per_female * s$max_females_per_male)
  }
})

test_that("suggested caps satisfy the feasibility floor", {
  for (SR in c(1, 2, 4)) for (MM in c(1, 2, 4)) {
    if (MM < (SR + 1) / 2) next
    s <- suggest_mate_caps(400, SR, MM)
    expect_silent(check_target_viability(
      sr_target = SR, mean_mates_target = MM,
      max_males_per_female = s$max_males_per_female,
      max_females_per_male = s$max_females_per_male))
  }
})

test_that("infeasible targets are refused with the boundary message", {
  expect_error(suggest_mate_caps(400, sr_target = 2, mean_mates_target = 1),
               "BOUNDARY ERROR")
  expect_error(suggest_mate_caps(400, sr_target = 4, mean_mates_target = 2),
               "at least 1.0 mate")
})

test_that("headroom below one is refused", {
  expect_error(suggest_mate_caps(400, 2, 2, headroom = 0.5), "below the mating success")
})

test_that("the single-configuration case is announced, not left to surprise", {
  # SR = 1, MM = 1 at no headroom admits only 1:1 blocks, so R-hat and ESS
  # come back undefined from an exact answer.
  expect_message(s <- suggest_mate_caps(200, 1, 1, headroom = 1), "single configuration")
  expect_equal(s$block_types, 1)
})

test_that("a thin iteration budget warns against the state space", {
  expect_warning(suggest_mate_caps(400, 1, 4, headroom = 4, n_iter = 1e5),
                 "per type")
  expect_silent(suggest_mate_caps(400, 2, 2, headroom = 2, n_iter = 1e6))
})

test_that("the suggested weight follows the linear rule", {
  expect_equal(suggest_mate_caps(100, 2, 2)$suggested_weight, 50)
  expect_equal(suggest_mate_caps(1600, 2, 2)$suggested_weight, 800)
})
