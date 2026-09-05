#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib mateR2, .registration = TRUE
## usethis namespace: end
NULL

# Global variables declaration to suppress R CMD check NOTEs caused by
# ggplot2 aesthetic mappings (aes)
utils::globalVariables(c(
  "Female", "Females", "MAP_Count", "Male", "Males", "Np",
  "Plot_Value", "Weight", "id", "iteration", "log_prob", "mm",
  "sr", "type", "x", "x_female", "x_male", "y", "y_female", "y_male"
))
