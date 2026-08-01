#' Summarize Mating, Sibling, and Parental Class Statistics Across a Pedigree Ensemble
#'
#' Converts a list of simulated pedigrees back to realized breeding matrices via \code{\link{ped2mat}}
#' and computes summary statistics across all posterior draws.
#'
#' @param ped_list A list of pedigree data frames returned by \code{\link{simulate_pedigree_ensemble}}.
#' @return A list containing three integrated data frames:
#'   \item{mating_stats}{Data frame of mating pattern metrics per draw (\code{\link{mat.stats}}).}
#'   \item{sibling_stats}{Data frame of sibling dyad and family metrics per draw (\code{\link{sib.stats}}).}
#'   \item{parental_classes}{List of parental yield class data frames per draw (\code{\link{parent.class.stats}}).}
#' @export
summarize_pedigree_ensemble <- function(ped_list) {
  if (!is.list(ped_list)) {
    stop("Input 'ped_list' must be a list of pedigree data frames.")
  }

  # Ensure clean element names
  if (is.null(names(ped_list)) || any(names(ped_list) == "")) {
    nms <- paste0("Draw_", seq_along(ped_list))
    if (!is.null(names(ped_list)) && names(ped_list)[1] == "MAP") {
      nms[1] <- "MAP"
    }
    names(ped_list) <- nms
  }

  mating_list  <- list()
  sibling_list <- list()
  parent_list  <- list()

  for (i in seq_along(ped_list)) {
    draw_name <- names(ped_list)[i]
    ped       <- ped_list[[i]]

    if (!is.data.frame(ped)) {
      stop(sprintf("Element '%s' in 'ped_list' is not a valid data frame.", draw_name))
    }

    # Step 1: Long Pedigree -> Realized Breeding Matrix
    realized_mat <- ped2mat(ped)

    # Step 2: Calculate Network & Mating Stats
    m_stat <- mat.stats(realized_mat)
    m_stat <- cbind(Draw = draw_name, m_stat)
    mating_list[[draw_name]] <- m_stat

    # Step 3: Calculate Sibling Dyads & Distributions
    s_stat <- sib.stats(realized_mat)
    s_stat <- cbind(Draw = draw_name, s_stat)
    sibling_list[[draw_name]] <- s_stat

    # Step 4: Calculate Parental Classes (Singletons, Doubletons, etc.)
    p_stat <- parent.class.stats(realized_mat)
    parent_list[[draw_name]] <- p_stat
  }

  # Combine into single tidy data frames
  mating_df  <- do.call(rbind, mating_list)
  sibling_df <- do.call(rbind, sibling_list)

  rownames(mating_df)  <- NULL
  rownames(sibling_df) <- NULL

  return(list(
    mating_stats     = mating_df,
    sibling_stats    = sibling_df,
    parental_classes = parent_list
  ))
}
