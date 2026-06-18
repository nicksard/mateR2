#' @title Calculate a Closeness Score
#' @description This function calculates a score representing how close an
#'   `actual` value is to a `target_value`. The score is defined as the minimum
#'   of the ratio (`actual / target_value`) and its inverse. This ensures the
#'   score is always between 0 and 1, with values closer to 1 indicating a
#'   better match. The function guarantees a non-zero, non-negative output for
#'   logarithmic stability in subsequent calculations.
#' @param actual The value being evaluated. Must be a positive, finite number.
#' @param target_value The target value for comparison. Must be a positive,
#'   finite number.
#' @return A numeric score (`0 <= score <= 1`) representing the closeness.
#'   Returns 0 if inputs are invalid.
#' @examples
#' # A perfect match
#' calculate_closeness_score(100, 100)
#' # A good match
#' calculate_closeness_score(90, 100)
#' # A poor match
#' calculate_closeness_score(10, 100)
calculate_closeness_score <- function(actual, target_value) {
  if (target_value <= 0 || !is.finite(target_value)) { return(0) }
  if (actual <= 0 || !is.finite(actual)) { return(0) }
  ratio <- actual / target_value
  inv_ratio <- target_value / actual
  score <- min(ratio, inv_ratio)
  return(max(score, .Machine$double.eps))
}
