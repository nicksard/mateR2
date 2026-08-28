#' @title Calculate Mating Pattern Statistics from a Breeding Matrix
#' @description This function analyzes a breeding matrix and calculates key
#'   statistics related to mating patterns, such as population size, mate counts,
#'   and reproductive success. It is a key tool for verifying the output of simulation methods.
#' @param mat The breeding matrix to be analyzed. It is expected to be a
#'   numeric matrix where rows represent males and columns represent females.
#'   Non-zero values indicate a mating event, and the value itself represents
#'   the number of offspring (Reproductive Success) from that mating.
#' @details
#' `mean.male.rs` and `mean.female.rs` are k-bar in the Crow & Kimura (1970) and
#' Caballero (1994) sense: offspring per parent computed within each sex. They
#' satisfy `mean.female.rs = SR * mean.male.rs` by construction. `var.male.rs`
#' and `var.female.rs` are the corresponding V_k, and the scaled forms are
#' V_k / k-bar, the quantity conventionally reported when comparing reproductive
#' skew against empirical systems.
#'
#' `offspring.per.parent` pools both sexes. It is NOT k-bar: every offspring is
#' counted once but divided over two parents, so it is half the per-sex mean when
#' the sex ratio is even. It was previously named `overall.mean.rs`.
#'
#' @return A data frame containing population sizes (`num.males`, `num.females`),
#'   the total count of mate pairs (`mp.count`), the mean, min and max values for
#'   mating success, and per-sex reproductive success summaries: mean (k-bar),
#'   variance (V_k) and scaled variance (V_k / k-bar).
#' @examples
#' # Assume 'my_matrix' is a breeding matrix (rows=males, cols=females)
#' # my_matrix <- matrix(c(1,2,0,0, 0,0,1,3), nrow=2, byrow=TRUE)
#' # mat_statistics <- mat.stats(my_matrix)
#' # print(mat_statistics)
#' @export
mat.stats <- function(mat) {

  # --- 0. Zero-Class Safety Filter ---
  # Count how many individuals failed to breed
  zero_males <- sum(rowSums(mat) == 0)
  zero_females <- sum(colSums(mat) == 0)

  # If any are found, warn the user before removing them
  if (zero_males > 0 || zero_females > 0) {
    warning(sprintf(
      "Dropped %d male(s) and %d female(s) with zero reproductive success. Network statistics are strictly calculated for the successful breeding pool (Np).",
      zero_males, zero_females
    ))
  }

  # Apply the filter to ensure we only evaluate successful breeders
  mat <- mat[rowSums(mat) > 0, colSums(mat) > 0, drop = FALSE]

  # --- 1. Population Sizes ---
  num.males <- nrow(mat)
  num.females <- ncol(mat)

  # --- 2. Mating Success (Unique Partners) ---
  # Create a binary matrix to count unique mates (ignoring offspring counts)
  mat_bin <- mat
  mat_bin[mat_bin > 0] <- 1

  mp.count <- sum(mat_bin) # Total number of breeding pairs

  male_mates <- rowSums(mat_bin)
  female_mates <- colSums(mat_bin)

  mean_m_mates <- mean(male_mates)
  mean_f_mates <- mean(female_mates)

  # --- 3. Reproductive Success (Total Offspring) ---
  # Use the raw matrix to sum the actual offspring values
  male_rs <- rowSums(mat)
  female_rs <- colSums(mat)
  total_offspring <- sum(mat)

  # --- 4. Compile Statistics ---
  stats_df <- data.frame(
    num.males = num.males,
    num.females = num.females,
    mp.count = mp.count,

    # Mating metrics
    mean.male.mates = round(mean_m_mates, 2),
    min.male.mates = min(male_mates),
    max.male.mates = max(male_mates),

    mean.female.mates = round(mean_f_mates, 2),
    min.female.mates = min(female_mates),
    max.female.mates = max(female_mates),

    # --- FIXED BIOLOGICAL MATH: Average of Averages ---
    overall.mean.mates = round((mean_m_mates + mean_f_mates) / 2.0, 2),

    # Reproductive Success (RS) metrics. These two are k-bar, per sex.
    mean.male.rs = round(mean(male_rs), 2),
    mean.female.rs = round(mean(female_rs), 2),

    # Variance in reproductive success (V_k) and its scaled form V_k / k-bar
    var.male.rs = round(stats::var(male_rs), 2),
    var.female.rs = round(stats::var(female_rs), 2),
    scaled.vk.male = round(stats::var(male_rs) / mean(male_rs), 2),
    scaled.vk.female = round(stats::var(female_rs) / mean(female_rs), 2),

    # Offspring per parent, pooled across both sexes. NOT k-bar: each offspring
    # is counted once but divided over two parents.
    offspring.per.parent = round(total_offspring / (num.males + num.females), 2)
  )

  return(stats_df)
}
