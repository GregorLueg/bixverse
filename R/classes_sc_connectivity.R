# single cell connectivity -----------------------------------------------------

# the neighbour graph lives in Rust behind an external pointer. the R object is
# a plain list around that pointer, so that the provenance stamp (an attribute)
# lands on the list and never on the externalptr itself. the pointer does not
# survive saveRDS() / qs2; save_sc_exp_to_disk() swaps it for bytes and
# load_existing() swaps it back.

## constructors ----------------------------------------------------------------

#' Wrap a Rust neighbour graph pointer
#'
#' @description
#' Low level constructor around the external pointer returned by the Rust
#' graph builders. Reads the shape back from Rust, so a graph built here and one
#' rehydrated from bytes end up identical.
#'
#' @param ptr External pointer to a Rust `SparseGraph<f32>`.
#' @param method String. How the graph was built, e.g. `"snn"`, `"umap"`,
#' `"knn"`, `"bbknn"`. `NA` if unknown.
#'
#' @returns A `SingleCellConnectivity` list with
#' \itemize{
#'   \item ptr - The external pointer.
#'   \item n - Integer. Number of nodes (cells).
#'   \item n_edges - Numeric. Number of undirected edges.
#'   \item method - String. How the graph was built.
#' }
#'
#' @keywords internal
.new_sc_connectivity <- function(ptr, method = NA_character_) {
  # checks
  checkmate::assertClass(ptr, "externalptr")
  checkmate::qassert(method, "s1")

  info <- rs_graph_info(ptr)
  if (info$directed) {
    stop("The neighbour graph must be undirected.")
  }

  res <- list(
    ptr = ptr,
    n = info$n_nodes,
    n_edges = info$n_edges,
    method = method
  )
  class(res) <- "SingleCellConnectivity"

  res
}

#' Convert a graph into a SingleCellConnectivity
#'
#' @description
#' Moves a neighbour graph into Rust. An igraph is flattened into an edge list
#' (edge attribute `weight` if present, otherwise all weights 1) and its graph
#' attribute `method` is carried over. A `SingleCellConnectivity` is returned
#' as is.
#'
#' @param x An undirected `igraph` or a `SingleCellConnectivity`.
#' @param ... Other parameters.
#'
#' @returns A `SingleCellConnectivity`.
#'
#' @export
#'
#' @examples
#' g <- igraph::make_graph("Zachary")
#' conn <- as_sc_connectivity(g)
#' conn
as_sc_connectivity <- function(x, ...) {
  UseMethod("as_sc_connectivity")
}

#' @rdname as_sc_connectivity
#'
#' @export
as_sc_connectivity.igraph <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "igraph")
  checkmate::assertFALSE(igraph::is_directed(x))

  edges <- igraph::as_edgelist(x, names = FALSE)
  weights <- igraph::E(x)$weight
  if (is.null(weights)) {
    weights <- rep(1, nrow(edges))
  }
  method <- igraph::graph_attr(x, "method")
  if (is.null(method)) {
    method <- NA_character_
  }

  ptr <- rs_graph_from_edges(
    edges = as.integer(t(edges)),
    weights = as.numeric(weights),
    n_nodes = igraph::vcount(x)
  )

  .new_sc_connectivity(ptr, method = as.character(method))
}

#' @rdname as_sc_connectivity
#'
#' @export
as_sc_connectivity.SingleCellConnectivity <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "SingleCellConnectivity")

  x
}

## primitives ------------------------------------------------------------------

#' @export
print.SingleCellConnectivity <- function(x, ...) {
  alive <- isTRUE(.connectivity_alive(x))
  cat(
    sprintf(
      "SingleCellConnectivity: %i cells, %.0f edges\n",
      x$n,
      x$n_edges
    ),
    sprintf("  Method: %s\n", x$method),
    if (!alive) "  Pointer: dead (reload the object)\n" else NULL,
    sep = ""
  )
  invisible(x)
}

#' Is the Rust pointer of a SingleCellConnectivity still valid?
#'
#' @param x A `SingleCellConnectivity`.
#'
#' @returns `TRUE`, or a string describing why the pointer cannot be used.
#'
#' @keywords internal
.connectivity_alive <- function(x) {
  ptr <- x[["ptr"]]
  if (is.null(ptr)) {
    return("The graph is in its on-disk form (bytes) and has no pointer.")
  }
  if (identical(ptr, methods::new("externalptr"))) {
    return(paste(
      "The graph pointer is dead, most likely because the object went",
      "through saveRDS()/readRDS() or qs2. Use save_sc_exp_to_disk() and",
      "load_existing(), or to_bytes() / from_bytes()."
    ))
  }
  TRUE
}

#' Get the live Rust pointer of a SingleCellConnectivity
#'
#' @param x A `SingleCellConnectivity`.
#'
#' @returns The external pointer. Errors if it is dead.
#'
#' @keywords internal
.connectivity_ptr <- function(x) {
  # checks
  checkmate::assertClass(x, "SingleCellConnectivity")

  alive <- .connectivity_alive(x)
  if (!isTRUE(alive)) {
    stop(alive)
  }

  x[["ptr"]]
}

## conversion ------------------------------------------------------------------

#' Convert a graph to igraph
#'
#' @description
#' Rebuilds the neighbour graph as an undirected igraph with the weights in
#' the edge attribute `weight` and the build method in the graph attribute
#' `method`. Costly on large graphs; the package's own consumers work on the
#' Rust graph directly.
#'
#' @param x A `SingleCellConnectivity` or an `igraph` (returned as is).
#' @param ... Other parameters.
#'
#' @returns An undirected, weighted `igraph`.
#'
#' @export
#'
#' @examples
#' conn <- as_sc_connectivity(igraph::make_graph("Zachary"))
#' as_igraph(conn)
as_igraph <- function(x, ...) {
  UseMethod("as_igraph")
}

#' @rdname as_igraph
#'
#' @export
as_igraph.SingleCellConnectivity <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "SingleCellConnectivity")

  edges <- rs_graph_to_edges(.connectivity_ptr(x))

  g <- igraph::make_empty_graph(n = x$n, directed = FALSE)
  g <- igraph::add_edges(g, edges$edges, attr = list(weight = edges$weights))
  if (!is.na(x$method)) {
    igraph::graph_attr(g, "method") <- x$method
  }

  g
}

#' @rdname as_igraph
#'
#' @export
as_igraph.igraph <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "igraph")

  x
}

#' Serialise an object to bytes
#'
#' @description
#' Serialises the Rust neighbour graph of a `SingleCellConnectivity` into a raw
#' vector, which survives `saveRDS()` and qs2 where the pointer does not. The
#' build method rides along as the attribute `method`. Read it back with
#' [from_bytes()].
#'
#' @param x A `SingleCellConnectivity`.
#' @param ... Other parameters.
#'
#' @returns A raw vector.
#'
#' @export
#'
#' @examples
#' conn <- as_sc_connectivity(igraph::make_graph("Zachary"))
#' bytes <- to_bytes(conn)
#' from_bytes(bytes)
to_bytes <- function(x, ...) {
  UseMethod("to_bytes")
}

#' @rdname to_bytes
#'
#' @export
to_bytes.SingleCellConnectivity <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "SingleCellConnectivity")

  bytes <- rs_graph_to_bytes(.connectivity_ptr(x))
  attr(bytes, "method") <- x$method

  bytes
}

#' Read a SingleCellConnectivity back from bytes
#'
#' @description
#' Inverse of [to_bytes()]. Bytes from a different structure, an older byte
#' layout or another float width error rather than decode into garbage.
#'
#' @param bytes Raw vector written by [to_bytes()].
#'
#' @returns A `SingleCellConnectivity`.
#'
#' @export
#'
#' @examples
#' conn <- as_sc_connectivity(igraph::make_graph("Zachary"))
#' from_bytes(to_bytes(conn))
from_bytes <- function(bytes) {
  # checks
  checkmate::assertRaw(bytes)

  method <- attr(bytes, "method", exact = TRUE)
  if (is.null(method)) {
    method <- NA_character_
  }

  .new_sc_connectivity(rs_graph_from_bytes(bytes), method = method)
}

## persistence -----------------------------------------------------------------

#' Swap the Rust pointer of a cached graph for bytes
#'
#' @description
#' Used by `save_sc_exp_to_disk()`. Rewrites the elements of the list only, so
#' the provenance stamp attribute is kept. A graph already in its on-disk form
#' is returned as is.
#'
#' @param x The cached graph, or `NULL`.
#'
#' @returns The graph with `bytes` in place of `ptr`.
#'
#' @keywords internal
.connectivity_to_disk <- function(x) {
  # checks
  checkmate::assertClass(x, "SingleCellConnectivity", null.ok = TRUE)

  if (is.null(x) || is.null(x[["ptr"]])) {
    return(x)
  }
  # `[[<-` keeps the attributes, so the stamp survives
  x[["bytes"]] <- rs_graph_to_bytes(.connectivity_ptr(x))
  x[["ptr"]] <- NULL

  x
}

#' Rehydrate a cached graph loaded off disk
#'
#' @description
#' Used by `load_existing()`, through the cache migrations. A graph in its
#' on-disk form gets its pointer back. A legacy igraph is converted into a
#' `SingleCellConnectivity` and dropped. The provenance stamp is kept in both
#' cases.
#'
#' @param x The cached graph, or `NULL`.
#'
#' @returns A `SingleCellConnectivity`, or `x` as is.
#'
#' @keywords internal
.connectivity_from_disk <- function(x) {
  # checks
  checkmate::assertMultiClass(
    x,
    c("SingleCellConnectivity", "igraph"),
    null.ok = TRUE
  )

  if (inherits(x, "igraph")) {
    stamp <- .read_stamp(x)
    res <- as_sc_connectivity(.drop_stamp(x))
    if (!is.null(stamp)) {
      res <- .set_stamp(res, stamp)
    }
    return(res)
  }
  if (inherits(x, "SingleCellConnectivity") && !is.null(x[["bytes"]])) {
    x[["ptr"]] <- rs_graph_from_bytes(x[["bytes"]])
    x[["bytes"]] <- NULL
  }

  x
}
