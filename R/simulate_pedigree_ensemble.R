#' Simulate Pedigrees Across a Posterior Ensemble of Count Vectors
#'
#' Takes a list of Mate-Pair Summary Table count vectors, expands each into an Individual-by-Individual
#' Matrix (IIM), applies Curveball degree-preserving rewiring, populates offspring fecundities,
#' and sub-samples juveniles to produce ready-to-use pedigrees.
#'
#' @param ensemble An object returned by \code{\link{sample_posterior_ensemble}} or \code{\link{generate_map_table}}.
#' @param mixing_I Numeric. Curveball rewiring intensity between 0.0 (block-diagonal) and 1.0 (full panmictic mixing). Default is 0.25.
#' @param min_fecundity Integer. Minimum offspring allocated per active pair-bond (default: 20).
#' @param max_fecundity Integer. Maximum offspring allocated per active pair-bond (default: 100).
#' @param fecundity_type Character. Fitness distribution type for \code{\link{brd.mat.fitness}} (default: \code{"uniform"}).
#' @param juvenile_sample_size Integer. Number of juveniles to sub-sample per pedigree (default: 500).
#' @param include_map Logical. Whether to include the MAP pedigree as the first element of the output list (default: \code{TRUE}).
#'
#' @return A list of ground-truth pedigree data frames with columns \code{c("off", "mom", "dad")}.
#' @export
simulate_pedigree_ensemble <- function(ensemble,
                                       mixing_I = 0.25,
                                       min_fecundity = 20,
                                       max_fecundity = 100,
                                       fecundity_type = "uniform",
                                       juvenile_sample_size = 500,
                                       include_map = TRUE) {

  # --- 1. Base Grid Extraction ---
  map_tab <- NULL
  if (!is.null(ensemble$map_c_k) && is.data.frame(ensemble$map_c_k)) {
    map_tab <- ensemble$map_c_k
  } else if (!is.null(ensemble$map_table) && is.data.frame(ensemble$map_table)) {
    map_tab <- ensemble$map_table
  } else {
    stop("Input 'ensemble' must contain a valid 'map_c_k' or 'map_table' data frame.")
  }

  if (!all(c("Males", "Females") %in% names(map_tab))) {
    stop("The map table in 'ensemble' must contain 'Males' and 'Females' columns.")
  }
  config_base <- map_tab[, c("Males", "Females")]

  # --- 2. Ensemble Vector List Extraction ---
  if (!is.null(ensemble$ensemble_c_k) && is.list(ensemble$ensemble_c_k)) {
    c_k_list <- ensemble$ensemble_c_k
  } else {
    stop("Input 'ensemble' does not contain a valid 'ensemble_c_k' list of sample vectors.")
  }

  # Ensure all input draw list elements have explicit names
  draw_names <- paste0("Draw_", seq_along(c_k_list))
  names(c_k_list) <- draw_names

  # Prepend MAP count vector if requested
  if (include_map) {
    map_vec <- NULL
    if ("MAP_Count" %in% names(map_tab)) {
      map_vec <- map_tab$MAP_Count
    } else if ("Count" %in% names(map_tab)) {
      map_vec <- map_tab$Count
    }

    if (!is.null(map_vec)) {
      c_k_list <- c(list(MAP = map_vec), c_k_list)
    }
  }

  # --- 3. Iterate over draws and simulate pedigrees ---
  pedigrees <- lapply(seq_along(c_k_list), function(i) {

    current_mp_table <- config_base
    current_mp_table$MAP_Count <- as.numeric(c_k_list[[i]])

    # Stage 1 -> Stage 2 Expansion (Mate-Pair Table to Binary Matrix)
    binary_mat <- mp_table_to_matrix(current_mp_table)

    # Stage 2 Curveball Rewiring (Passing parameter I = mixing_I)
    if (mixing_I > 0) {
      mixed_mat <- randomize_mating_structure(binary_mat, I = mixing_I)
    } else {
      mixed_mat <- binary_mat
    }

    # Stage 3 Fecundity Assignment
    fecund_mat <- brd.mat.fitness(
      mat      = mixed_mat,
      min.fert = min_fecundity,
      max.fert = max_fecundity,
      type     = fecundity_type
    )

    # Sub-sample Juveniles
    sampled_df <- mat.sub.sample(fecund_mat, num_offspring = juvenile_sample_size)

    # Convert to standard Pedigree Data Frame
    ped <- convert2ped(sampled_df)
    return(ped)
  })

  # Assign clean list names
  names(pedigrees) <- names(c_k_list)

  return(pedigrees)
}
