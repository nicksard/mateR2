#' @title Count Connected Components of a Bipartite Mating Network
#'
#' @description Quantifies macro-topological fragmentation by converting a binary
#'   mating matrix into an undirected bipartite graph and counting its connected
#'   components. Matrices expanded straight from a Mate-Pair Summary Table are
#'   block-diagonal and fragment into many isolated components; as
#'   \code{\link{randomize_mating_structure}} rewires the network these dissolve
#'   toward a single integrated population graph.
#'
#' @details
#' Degree-preserving edge swapping holds every node degree constant, so
#' first-order properties such as connectance are invariant under rewiring.
#' Component counts capture the higher-order reorganisation that rewiring does
#' produce. Isolated nodes each count as their own component, matching igraph's
#' convention.
#'
#' Requires the igraph package, which is listed under Suggests.
#'
#' @param mat A binary mating matrix (rows = males, columns = females). Non-zero
#'   entries are treated as edges, so a weighted breeding matrix may also be
#'   passed and will be binarised.
#'
#' @return An integer: the number of connected components.
#'
#' @seealso \code{\link{comating_pair_density}}
#'
#' @examples
#' \dontrun{
#' block <- mp_table_to_matrix(data.frame(Males = 2, Females = 2, Count = 25))
#' count_network_components(block)                                   # many
#' count_network_components(randomize_mating_structure(block, 0.25)) # few
#' }
#'
#' @export
count_network_components <- function(mat) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for count_network_components(). ",
         "Install it with install.packages('igraph').", call. = FALSE)
  }
  mat <- as.matrix(mat)
  mat_bin <- matrix(as.integer(mat != 0), nrow = nrow(mat), ncol = ncol(mat))

  # igraph renamed graph_from_incidence_matrix() to graph_from_biadjacency_matrix()
  # in 2.0.0; support both so this works across installed versions.
  if (exists("graph_from_biadjacency_matrix", where = asNamespace("igraph"))) {
    g <- igraph::graph_from_biadjacency_matrix(mat_bin)
  } else {
    g <- igraph::graph_from_incidence_matrix(mat_bin)
  }
  as.integer(igraph::count_components(g))
}
