#' @title Convert Offspring Data Frame to Pedigree Format
#' @description This function takes a data frame of sampled offspring (as output
#'   by `mat.sub.sample()`) and converts it into a long-format pedigree file,
#'   which is a standard input for genetic and parentage analysis software.
#' @param df A data frame with at least the columns `dads`, `moms`, and `off1`,
#'   representing father, mother, and the number of sampled offspring per pair.
#' @return A data frame with three columns: `off` (offspring ID), `mom` (mother ID),
#'   and `dad` (father ID), with one row for each sampled offspring.
#' @export
convert2ped <- function(df) {
  # Making sure there are no zeros
  df <- df[df$off1 != 0, ]
  
  # Check if the data frame is empty after filtering
  if (nrow(df) == 0) {
    warning("Input data frame contains no offspring. Returning an empty data frame.")
    return(data.frame(off = character(0), mom = character(0), dad = character(0)))
  }
  
  # Create a generic pedigree with the remaining offspring
  # Use an alternative to the loop for efficiency
  offspring_counts <- df$off1
  moms_vec <- df$moms
  dads_vec <- df$dads
  
  df.out <- data.frame(
    mom = rep(moms_vec, times = offspring_counts),
    dad = rep(dads_vec, times = offspring_counts)
  )
  
  # Add unique offspring IDs
  df.out$off <- paste0("off_", 1:nrow(df.out))
  
  # Reorder columns to the standard pedigree format
  df.out <- df.out[, c("off", "mom", "dad")]
  
  return(df.out)
}
