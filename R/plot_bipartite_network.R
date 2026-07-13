#' @title Plot a Bipartite Mating Network
#' @description Visualizes a breeding matrix as a bipartite network using ggplot2.
#' @param mat A numeric matrix representing the breeding population. Rows = males, Cols = females.
#' @param title A character string for the plot title.
#' @param is_weighted Logical. If TRUE, edge thickness scales with the matrix values (fecundity).
#' @return A `ggplot` object representing the bipartite network.
#' @import ggplot2
#' @export
plot_bipartite_network <- function(mat, title = "Bipartite Mating Network", is_weighted = FALSE) {

  # --- Input Validation and Setup ---
  mat <- as.matrix(mat)
  if (is.null(rownames(mat))) rownames(mat) <- paste0("M", 1:nrow(mat))
  if (is.null(colnames(mat))) colnames(mat) <- paste0("F", 1:ncol(mat))

  Nm <- nrow(mat)
  Nf <- ncol(mat)
  max_nodes <- max(Nm, Nf)

  # --- 1. Generate Node Coordinates ---
  # Males on x=1, Females on x=2.
  # Y-coordinates are reversed so M1/F1 are at the top, and centered if Nm != Nf.
  nodes_males <- data.frame(
    id = rownames(mat),
    x = 1,
    y = seq(Nm, 1) + (max_nodes - Nm) / 2,
    type = "Male"
  )
  nodes_females <- data.frame(
    id = colnames(mat),
    x = 2,
    y = seq(Nf, 1) + (max_nodes - Nf) / 2,
    type = "Female"
  )
  nodes <- rbind(nodes_males, nodes_females)

  # --- 2. Generate Edge Coordinates ---
  # Convert matrix to a long dataframe of edges
  mat_df <- as.data.frame(as.table(mat))
  colnames(mat_df) <- c("Male", "Female", "Weight")
  edges <- mat_df[mat_df$Weight > 0, ] # Only keep actual matings

  # Merge coordinates to find start (x,y) and end (x,y) for each edge
  edges <- merge(edges, nodes_males[, c("id", "x", "y")], by.x = "Male", by.y = "id")
  colnames(edges)[4:5] <- c("x_male", "y_male")
  edges <- merge(edges, nodes_females[, c("id", "x", "y")], by.x = "Female", by.y = "id")
  colnames(edges)[6:7] <- c("x_female", "y_female")

  # --- 3. Build the Plot ---
  p <- ggplot2::ggplot()

  # Draw Edges (Binary vs Weighted)
  if (is_weighted) {
    p <- p + ggplot2::geom_segment(data = edges, ggplot2::aes(x = x_male, y = y_male, xend = x_female, yend = y_female, linewidth = Weight), color = "gray60", alpha = 0.7) +
      ggplot2::scale_linewidth_continuous(range = c(0.5, 3), guide = "none") # Thicker lines for more offspring
  } else {
    p <- p + ggplot2::geom_segment(data = edges, ggplot2::aes(x = x_male, y = y_male, xend = x_female, yend = y_female), linewidth = 1, color = "gray60", alpha = 0.7)
  }

  # Draw Nodes and Labels
  p <- p +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x = x, y = y, color = type, shape = type), size = 6) +
    ggplot2::scale_color_manual(values = c("Male" = "#377eb8", "Female" = "#e41a1c")) +
    ggplot2::scale_shape_manual(values = c("Male" = 15, "Female" = 16)) + # 15=Square, 16=Circle
    # Nudge labels outside the graph
    ggplot2::geom_text(data = nodes, ggplot2::aes(x = ifelse(x == 1, x - 0.2, x + 0.2), y = y, label = id), size = 4, fontface = "bold") +

    # Styling
    ggplot2::xlim(0.5, 2.5) +
    ggplot2::theme_void() +
    ggplot2::labs(title = title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "none" # Hide legend for clean look
    )

  return(p)
}
