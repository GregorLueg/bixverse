# sc connectivity (Rust neighbour graph) ---------------------------------------

source("helper_sc.R", local = TRUE)

set.seed(123L)

test_temp_dir <- sc_test_dir("sc_connectivity")

## fixture ---------------------------------------------------------------------

single_cell_test_data <- sc_test_fixture()

sc_object <- sc_test_object(
  test_temp_dir,
  single_cell_test_data,
  sc_qc_param = sc_test_qc_params(single_cell_test_data, target_size = 1000)
)

sc_object <- find_hvg_sc(sc_object, hvg_no = 30L, .verbose = FALSE)
sc_object <- calculate_pca_sc(sc_object, no_pcs = 10L, .verbose = FALSE)
sc_object <- find_neighbours_sc(sc_object, .verbose = FALSE)

conn <- get_connectivities(sc_object)
n_cells <- length(get_cells_to_keep(sc_object))

# edges as a sorted (from, to, weight) table; edge order is not part of the
# contract
edge_table <- function(edges, weights) {
  pairs <- matrix(edges, ncol = 2L, byrow = TRUE)
  dt <- data.table::data.table(
    from = pmin(pairs[, 1], pairs[, 2]),
    to = pmax(pairs[, 1], pairs[, 2]),
    weight = weights
  )
  data.table::setorder(dt, from, to)
  dt
}

conn_edges <- function(x) {
  with(rs_graph_to_edges(bixverse:::.connectivity_ptr(x)), {
    edge_table(edges, weights)
  })
}

igraph_edges <- function(g) {
  edge_table(
    as.integer(t(igraph::as_edgelist(g, names = FALSE))),
    igraph::E(g)$weight
  )
}

## class -----------------------------------------------------------------------

expect_inherits(
  current = conn,
  class = "SingleCellConnectivity",
  info = "connectivity - find_neighbours_sc() stores the Rust graph"
)

expect_equal(
  current = conn$n,
  target = n_cells,
  info = "connectivity - one node per kept cell"
)

expect_equal(
  current = conn$method,
  target = "snn",
  info = "connectivity - build method recorded"
)

expect_stdout(
  current = print(conn),
  pattern = "SingleCellConnectivity",
  info = "connectivity - print method dispatches"
)

## igraph <-> ptr --------------------------------------------------------------

g <- as_igraph(conn)

expect_equal(
  current = igraph::vcount(g),
  target = n_cells,
  info = "as_igraph - vertex count"
)

expect_equal(
  current = igraph::ecount(g),
  target = conn$n_edges,
  info = "as_igraph - edge count matches the Rust graph"
)

expect_equal(
  current = igraph_edges(g),
  target = conn_edges(conn),
  info = "as_igraph - edge set identical to the Rust graph"
)

conn_back <- as_sc_connectivity(g)

expect_equal(
  current = conn_edges(conn_back),
  target = conn_edges(conn),
  info = "as_sc_connectivity - igraph -> ptr gives the same edge set"
)

expect_equal(
  current = conn_back$method,
  target = "snn",
  info = "as_sc_connectivity - method carried over from the graph attribute"
)

expect_identical(
  current = as_sc_connectivity(conn),
  target = conn,
  info = "as_sc_connectivity - identity on SingleCellConnectivity"
)

## leiden / louvain ------------------------------------------------------------

for (algo in c("leiden", "louvain")) {
  expect_identical(
    current = leiden_clustering(conn, algorithm = algo),
    target = leiden_clustering(g, algorithm = algo),
    info = sprintf("%s - identical labels via the pointer and via igraph", algo)
  )
}

sc_igraph <- set_connectivities(sc_object, g, from = "knn")

expect_inherits(
  current = get_connectivities(sc_igraph),
  class = "SingleCellConnectivity",
  info = "set_connectivities - igraph converted on entry"
)

expect_identical(
  current = find_clusters_sc(sc_igraph)[["leiden_clustering"]],
  target = find_clusters_sc(sc_object)[["leiden_clustering"]],
  info = "find_clusters_sc - same labels from an igraph set by hand"
)

## bytes -----------------------------------------------------------------------

bytes <- to_bytes(conn)

expect_true(
  current = is.raw(bytes),
  info = "to_bytes - raw vector"
)

conn_bytes <- from_bytes(bytes)

expect_equal(
  current = conn_edges(conn_bytes),
  target = conn_edges(conn),
  info = "bytes - round trip keeps the edge set"
)

expect_equal(
  current = conn_bytes$method,
  target = conn$method,
  info = "bytes - round trip keeps the method"
)

bad_bytes <- bytes
bad_bytes[1] <- as.raw(0)

expect_error(
  current = from_bytes(bad_bytes),
  pattern = "magic",
  info = "bytes - foreign bytes are rejected"
)

bad_version <- bytes
bad_version[9] <- as.raw(99)

expect_error(
  current = from_bytes(bad_version),
  pattern = "version",
  info = "bytes - another format version is rejected"
)

## dead pointer ----------------------------------------------------------------

conn_dead <- unserialize(serialize(conn, NULL))

expect_error(
  current = leiden_clustering(conn_dead),
  pattern = "pointer is dead",
  info = "dead pointer - clear error after serialisation"
)

## persistence -----------------------------------------------------------------

sc_object <- find_clusters_sc(sc_object)
labels <- sc_object[["leiden_clustering"]]

for (type in c("qs2", "rds")) {
  unlink(file.path(test_temp_dir, c("memory.qs2", "memory.rds")))
  save_sc_exp_to_disk(sc_object, type = type)

  # the object in memory keeps its live pointer
  expect_true(
    current = isTRUE(bixverse:::.connectivity_alive(
      get_connectivities(sc_object)
    )),
    info = sprintf("save (%s) - the live object is untouched", type)
  )

  reloaded <- load_existing(
    SingleCells(dir_data = test_temp_dir),
    .verbose = FALSE
  )
  reloaded_conn <- get_connectivities(reloaded)

  expect_inherits(
    current = reloaded_conn,
    class = "SingleCellConnectivity",
    info = sprintf(
      "load (%s) - graph comes back as SingleCellConnectivity",
      type
    )
  )

  expect_equal(
    current = conn_edges(reloaded_conn),
    target = conn_edges(conn),
    info = sprintf("load (%s) - edge set survives", type)
  )

  expect_false(
    current = any(get_sc_cache_status(reloaded)$stale),
    info = sprintf("load (%s) - provenance stamp survives", type)
  )

  expect_identical(
    current = find_clusters_sc(reloaded)[["leiden_clustering"]],
    target = labels,
    info = sprintf("load (%s) - clustering identical after reload", type)
  )
}

## legacy igraph cache ---------------------------------------------------------

unlink(file.path(test_temp_dir, c("memory.qs2", "memory.rds")))
save_sc_exp_to_disk(sc_object, type = "qs2")

memory_path <- file.path(test_temp_dir, "memory.qs2")
saved <- qs2::qs_read(memory_path)
stamp <- attr(saved$sc_cache$connectivities, "bixverse_stamp")
legacy_graph <- g
attr(legacy_graph, "bixverse_stamp") <- stamp
saved$sc_cache$connectivities <- legacy_graph
qs2::qs_save(saved, memory_path)

legacy <- load_existing(SingleCells(dir_data = test_temp_dir), .verbose = FALSE)
legacy_conn <- get_sc_cache(legacy)$connectivities

expect_inherits(
  current = legacy_conn,
  class = "SingleCellConnectivity",
  info = "legacy - an igraph cache is converted on load"
)

expect_equal(
  current = attr(legacy_conn, "bixverse_stamp"),
  target = stamp,
  info = "legacy - the stamp is carried over from the igraph"
)

expect_equal(
  current = conn_edges(legacy_conn),
  target = conn_edges(conn),
  info = "legacy - edge set identical after conversion"
)

expect_identical(
  current = find_clusters_sc(legacy)[["leiden_clustering"]],
  target = labels,
  info = "legacy - clustering identical after conversion"
)

## cleanup ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir)
