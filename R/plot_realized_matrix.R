#' @title Plot a Realized Breeding Matrix
#' @description Visualizes a specific breeding matrix (Individual-by-Individual) as a heatmap.
#'   Males are represented on the rows (Y-axis) and Females on the columns (X-axis).
#'   Zero values (no mating) are colored white, while successful matings are
#'   colored using a custom heat gradient.
#' @param mat A numeric matrix representing the breeding population. Rows = males, Cols = females.
#' @param title A character string for the plot title. Default is "Realized Breeding Matrix".
#' @return A `ggplot` object representing the matrix heatmap.
#' @import ggplot2
#' @export
#' @examples
#' # Create a sample 10x10 matrix with sparse matings
#' # test_mat <- matrix(sample(c(0, 0, 0, 1, 2, 3), 100, replace=TRUE), nrow=10)
#' # plot_realized_matrix(test_mat)
plot_realized_matrix <- function(mat, title = "Realized Breeding Matrix",fill_label = "Mates") {

  # --- Input Validation ---
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("Input 'mat' must be a matrix or data frame.")
  }

  # Ensure it is a matrix for dimnames handling
  mat <- as.matrix(mat)

  # --- Assign Dimension Names (if missing) ---
  if (is.null(rownames(mat))) rownames(mat) <- paste0("M", 1:nrow(mat))
  if (is.null(colnames(mat))) colnames(mat) <- paste0("F", 1:ncol(mat))

  # --- Data Transformation ---
  # Use base R to neatly melt the matrix into a long-format dataframe for ggplot
  mat_df <- as.data.frame(as.table(mat))
  colnames(mat_df) <- c("Male", "Female", "Offspring")

  # Convert 0s to NA so they can be mapped to purely white.
  # The gradient will only apply to cells with >0 offspring.
  mat_df$Plot_Value <- mat_df$Offspring
  mat_df$Plot_Value[mat_df$Plot_Value == 0] <- NA

  # Reverse the Male factor levels so Row 1 (M1) plots at the TOP of the graph
  mat_df$Male <- factor(mat_df$Male, levels = rev(rownames(mat)))

  # --- Generate Plot ---
  p <- ggplot2::ggplot(mat_df, ggplot2::aes(x = Female, y = Male, fill = Plot_Value)) +
    # Draw the tiles. A faint grey border helps separate the white cells slightly.
    ggplot2::geom_tile(color = "grey95", linewidth = 0.1) +

    # Apply your custom color palette for >0 values, and set NA (0s) to white
    ggplot2::scale_fill_gradientn(
      colors = c("#ffeda0", "#feb24c", "#f03b20"),
      na.value = "white",
      name = fill_label
    ) +

    # Styling
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = title,
      x = "Females (Columns)",
      y = "Males (Rows)",
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      # Shrink and rotate axis text in case of large matrices (e.g. 100x100)
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                          hjust=1, size = 6, color = "black"),
      axis.text.y = ggplot2::element_text(size = 6, color = "black"),

      # Put the legend on the right for standard heatmap feel
      legend.position = "right"
    )

  return(p)
}
