#' @title Plot Mate Pair Table Block Distribution
#' @description Visualizes a mate pair summary table (e.g., a MAP estimate from the MCMC sampler)
#'   as a heatmap. It maps Males to the Y-axis, Females to the X-axis, and colors the tiles
#'   based on the frequency of each mating configuration.
#' @param mp_table A data frame containing the mate pair counts. It must contain
#'   `Males`, `Females`, and `MAP_Count` columns.
#' @param title A character string for the plot title. Default is "Mating Matrix Block Distribution".
#' @return A `ggplot` object representing the heatmap.
#' @import ggplot2
#' @export
#' @examples
#' # sample_data <- data.frame(Males = c(1, 1, 2), Females = c(1, 2, 1), MAP_Count = c(16, 12, 0))
#' # plot_mate_matrix(sample_data)
plot_mate_matrix <- function(mp_table, title = "Mating Matrix Block Distribution") {

  # --- Input Validation ---
  if (!all(c("Males", "Females", "MAP_Count") %in% names(mp_table))) {
    stop("Input table must contain 'Males', 'Females', and 'MAP_Count' columns.")
  }

  # Safely duplicate the data frame
  plot_data <- as.data.frame(mp_table)

  # Convert to factors for crisp categorical plotting
  plot_data$Males <- factor(plot_data$Males)
  plot_data$Females <- factor(plot_data$Females)

  # Determine the midpoint for text color contrast
  max_val <- max(plot_data$MAP_Count, na.rm = TRUE)

  # --- Generate the Plot ---
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Females, y = Males, fill = MAP_Count)) +
    ggplot2::geom_tile(color = "gray50", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = MAP_Count),
      color = ifelse(plot_data$MAP_Count > (max_val / 2), "white", "black"),
      fontface = "bold", size = 5
    ) +
    # Add "white" to the start of the gradient so 0 counts are blank but keep their numbers!
    ggplot2::scale_fill_gradientn(
      colors = c("white", "#ffeda0", "#feb24c", "#f03b20"),
      name = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = title,
      x = "Max Females per Male (Columns)",
      y = "Max Males per Female (Rows)"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(face = "bold", color="black")
    )

  return(p)
}
