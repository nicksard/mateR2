#' @title Convert a Breeding Matrix to a Mate Pair Summary Table
#'
#' @description Recovers the Mate-Pair Summary Table from a binary bipartite
#'   mating matrix, inverting \code{\link{mp_table_to_matrix}}. Each mating block
#'   is a connected component of the bipartite graph, and the table reports how
#'   many blocks of each \emph{m}:\emph{f} type the matrix contains.
#'
#' @details
#' \strong{A block is a connected component, not a degree signature.} Two
#' individuals belong to the same mating block when a path of realized pairings
#' runs between them. The number of mates an individual has does not identify its
#' block: two separate 1:1 pairs and one 2:2 block contain individuals whose
#' degrees alone cannot tell them apart, and only the connectivity distinguishes
#' them. Components are found with a union-find pass over the edge list, which
#' needs no external graph package.
#'
#' \strong{This inverts Stage 1, not Stage 2.} \code{mp_table_to_matrix()}
#' expands a count vector into complete, block-diagonal groups, and this function
#' recovers that vector exactly. Once \code{\link{randomize_mating_structure}}
#' has rewired the matrix at any mixing intensity above zero, blocks merge and
#' the components stop being complete bipartite groups. The table then describes
#' the components actually present, which is a real property of the mixed network
#' but is not the block vector that generated it. A warning fires when incomplete
#' components are found, because the distinction is easy to lose.
#'
#' @param mat A binary matrix, rows males and columns females.
#' @param drop_unmated Logical. Individuals with no mates form components with no
#'   edges and are not valid mating blocks. \code{TRUE} (default) removes them
#'   and reports how many; \code{FALSE} keeps them as 1:0 and 0:1 rows.
#'
#' @return A data frame with columns \code{Males}, \code{Females} and
#'   \code{Count}, one row per distinct block type. Carries attributes
#'   \code{n_unmated_males}, \code{n_unmated_females} and \code{complete_blocks},
#'   the last being \code{FALSE} when any component is not a complete bipartite
#'   group.
#'
#' @seealso \code{\link{mp_table_to_matrix}} for the forward direction,
#'   \code{\link{count_network_components}} for the component count alone.
#'
#' @examples
#' tbl <- data.frame(Males = c(1, 2), Females = c(1, 2), Count = c(2, 1))
#' m   <- mp_table_to_matrix(tbl)
#' mat2mp_sum_table(m)   # recovers tbl
#'
#' @export
mat2mp_sum_table <- function(mat, drop_unmated = TRUE) {

  if (is.data.frame(mat)) mat <- as.matrix(mat)
  if (!is.matrix(mat) || !is.numeric(mat) || anyNA(mat) ||
      any(mat != 0 & mat != 1)) {
    stop("\n[mateR2 INPUT ERROR]\n'mat' must be a complete numeric matrix of 0s and 1s, rows males and columns females.",
         call. = FALSE)
  }

  n_m <- nrow(mat)
  n_f <- ncol(mat)
  if (n_m == 0L || n_f == 0L) {
    stop("\n[mateR2 INPUT ERROR]\n'mat' has a zero dimension; there is no mating structure to summarise.",
         call. = FALSE)
  }

  edges <- which(mat == 1, arr.ind = TRUE)
  if (nrow(edges) == 0L) {
    warning("'mat' contains no mating links, so there are no blocks to report.",
            call. = FALSE)
    out <- data.frame(Males = integer(0), Females = integer(0),
                      Count = integer(0))
    attr(out, "n_unmated_males")   <- n_m
    attr(out, "n_unmated_females") <- n_f
    attr(out, "complete_blocks")   <- TRUE
    return(out)
  }

  # --- Union-find over males 1..n_m and females n_m+1..n_m+n_f -------------
  parent <- seq_len(n_m + n_f)
  find <- function(a) {
    root <- a
    while (parent[root] != root) root <- parent[root]
    while (parent[a] != root) { nxt <- parent[a]; parent[a] <<- root; a <- nxt }
    root
  }
  for (e in seq_len(nrow(edges))) {
    ra <- find(edges[e, 1])
    rb <- find(n_m + edges[e, 2])
    if (ra != rb) parent[ra] <- rb
  }
  comp <- vapply(seq_len(n_m + n_f), find, integer(1))

  # --- Count males and females in each component ---------------------------
  ids     <- sort(unique(comp))
  males   <- as.integer(table(factor(comp[seq_len(n_m)],       levels = ids)))
  females <- as.integer(table(factor(comp[n_m + seq_len(n_f)], levels = ids)))

  # A block built by mp_table_to_matrix() is fully connected, so its edge count
  # is m * f. Rewiring breaks that, and the result is then a component summary
  # rather than the generating block vector.
  obs_edges <- as.integer(table(factor(comp[edges[, 1]], levels = ids)))
  complete  <- all(obs_edges == males * females)

  unmated_m <- sum(males == 1L & females == 0L)
  unmated_f <- sum(males == 0L & females == 1L)
  if (drop_unmated) {
    keep    <- males > 0L & females > 0L
    males   <- males[keep]
    females <- females[keep]
  }

  tab   <- table(paste(males, females, sep = ":"))
  parts <- do.call(rbind, strsplit(names(tab), ":", fixed = TRUE))
  out <- data.frame(Males   = as.integer(parts[, 1]),
                    Females = as.integer(parts[, 2]),
                    Count   = as.integer(tab),
                    stringsAsFactors = FALSE)
  out <- out[order(out$Males, out$Females), ]
  rownames(out) <- NULL

  if (!complete) {
    warning("Some components are not complete bipartite blocks, which is what rewiring produces. The table describes the components present, not the block vector that generated the matrix.",
            call. = FALSE)
  }
  if (drop_unmated && (unmated_m + unmated_f) > 0L) {
    message(sprintf("Dropped %d unmated male(s) and %d unmated female(s); they form no mating block.",
                    unmated_m, unmated_f))
  }
  attr(out, "n_unmated_males")   <- unmated_m
  attr(out, "n_unmated_females") <- unmated_f
  attr(out, "complete_blocks")   <- complete
  out
}
