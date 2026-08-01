# spatial graph construction ---------------------------------------------------

# Complements test_sp_analysis.R, which exercises the hex lattice through a
# SpatialSpotSubset. Here the target is `rs_sp_build_graph()` itself: every
# layout, every weighting, checked against a plain R reference, plus how
# `build_spatial_graph_sp()` fans the thing out over several samples.

source("helper_sp.R", local = TRUE)

test_temp_dir <- file.path(tempdir(), "sp_graph")
dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

## fixtures --------------------------------------------------------------------

n_rows <- 8L
n_cols <- 8L

hex_grid <- sp_lattice_positions("hex", n_rows, n_cols)
square_grid <- sp_lattice_positions("square", n_rows, n_cols)

as_coords <- function(grid) {
  out <- matrix(
    c(
      as.numeric(grid$pxl_col_in_fullres),
      as.numeric(grid$pxl_row_in_fullres)
    ),
    ncol = 2L
  )
  colnames(out) <- c("x", "y")
  out
}

hex_coords <- as_coords(hex_grid)
square_coords <- as_coords(square_grid)
n_spots <- nrow(hex_coords)

# `rs_sp_build_graph()` takes the array coordinates inside the parameter list.
# `params_sp_graph()` deliberately does not carry them, the method injects them
# from the obs table, so the tests do the same by hand.
with_array <- function(graph_params, grid) {
  graph_params$array_row <- as.integer(grid$array_row)
  graph_params$array_col <- as.integer(grid$array_col)
  graph_params
}

build <- function(coords, graph_params, grid = NULL, knn_params = NULL) {
  if (!is.null(grid)) {
    graph_params <- with_array(graph_params, grid)
  }
  rs_sp_build_graph(
    coords = coords,
    graph_params = graph_params,
    knn_params = knn_params,
    seed = 42L,
    verbose = 0L
  )
}

# pairwise distances, for the weighting references
hex_dist <- as.matrix(stats::dist(hex_coords))
diag(hex_dist) <- Inf

## the returned shape ----------------------------------------------------------

hex_binary <- build(hex_coords, params_sp_graph("hex", "binary"), hex_grid)

expect_equal(
  current = names(hex_binary),
  target = c("indices", "weights", "degree"),
  info = "sp graph - rs_sp_build_graph returns indices, weights and degree"
)

expect_equal(
  current = c(
    length(hex_binary$indices),
    length(hex_binary$weights),
    length(hex_binary$degree)
  ),
  target = rep(n_spots, 3L),
  info = "sp graph - one entry per spot in all three slots"
)

expect_equal(
  current = lengths(hex_binary$weights),
  target = lengths(hex_binary$indices),
  info = "sp graph - indices and weights are aligned element for element"
)

expect_true(
  current = all(vapply(
    hex_binary$indices,
    \(idx) !is.unsorted(idx),
    logical(1)
  )),
  info = "sp graph - neighbour lists come back sorted ascending"
)

expect_true(
  current = max(unlist(hex_binary$indices)) == n_spots - 1L &&
    min(unlist(hex_binary$indices)) == 0L,
  info = "sp graph - neighbour indices are local and 0-based"
)

## degree is not the neighbour count -------------------------------------------

# The single easiest thing to get wrong here. `degree` is the row sum of the
# symmetrised graph W + t(W), not the number of neighbours. On a reciprocal
# binary lattice that is exactly twice the neighbour count, and anywhere else
# it is neither.

hex_nb_counts <- lengths(hex_binary$indices)

expect_equal(
  current = as.numeric(hex_binary$degree),
  target = as.numeric(2L * hex_nb_counts),
  info = "sp graph - degree is twice the neighbour count on a binary lattice"
)

expect_equal(
  current = stats::median(hex_nb_counts),
  target = 6,
  info = "sp graph - the median hex spot has six neighbours"
)

expect_equal(
  current = stats::median(hex_binary$degree),
  target = 12,
  info = "sp graph - and a median symmetrised degree of twelve"
)

# the same identity through the sparse matrix the class actually caches
hex_sparse <- bixverse:::.sp_adjacency_to_sparse(
  hex_binary$indices,
  hex_binary$weights,
  n_spots
)

expect_equal(
  current = as.numeric(Matrix::rowSums(hex_sparse + Matrix::t(hex_sparse))),
  target = as.numeric(hex_binary$degree),
  info = "sp graph - degree equals the row sums of W + t(W)"
)

## lattice layouts against an R reference --------------------------------------

# The lattice layouts are supposed to be an exact table look-up over the array
# coordinates. `sp_lattice_neighbours()` does that look-up in plain R, so it is
# an independent reference rather than a restatement of the Rust.

expect_equal(
  current = lapply(hex_binary$indices, as.integer),
  target = sp_lattice_neighbours(
    hex_grid$array_row,
    hex_grid$array_col,
    "hex"
  ),
  info = "sp graph - the hex adjacency matches the honeycomb offsets"
)

square_4 <- build(
  square_coords,
  params_sp_graph("square", "binary", connectivity = 4L),
  square_grid
)

expect_equal(
  current = lapply(square_4$indices, as.integer),
  target = sp_lattice_neighbours(
    square_grid$array_row,
    square_grid$array_col,
    "square",
    4L
  ),
  info = "sp graph - the 4-connected square adjacency matches von Neumann"
)

expect_equal(
  current = max(lengths(square_4$indices)),
  target = 4L,
  info = "sp graph - a 4-connected interior bin has four neighbours"
)

square_8 <- build(
  square_coords,
  params_sp_graph("square", "binary", connectivity = 8L),
  square_grid
)

expect_equal(
  current = lapply(square_8$indices, as.integer),
  target = sp_lattice_neighbours(
    square_grid$array_row,
    square_grid$array_col,
    "square",
    8L
  ),
  info = "sp graph - the 8-connected square adjacency matches Moore"
)

expect_equal(
  current = max(lengths(square_8$indices)),
  target = 8L,
  info = "sp graph - an 8-connected interior bin has eight neighbours"
)

## knn and radius --------------------------------------------------------------

knn_graph <- build(hex_coords, params_sp_graph("knn", "binary", k = 5L))

expect_equal(
  current = unique(lengths(knn_graph$indices)),
  target = 5L,
  info = "sp graph - the knn layout gives every spot exactly k neighbours"
)

expect_true(
  current = !all(
    vapply(
      seq_len(n_spots),
      \(i) {
        setequal(knn_graph$indices[[i]], hex_binary$indices[[i]])
      },
      logical(1)
    )
  ),
  info = "sp graph - a k = 5 knn is not the same graph as the hex lattice"
)

# in-row hex neighbours sit 100 px apart, the out-of-row ones 100.34 px, so a
# 110 px cut-off takes the full honeycomb and nothing beyond it
radius_graph <- build(
  hex_coords,
  params_sp_graph("radius", "binary", radius = 110)
)

expect_equal(
  current = lapply(radius_graph$indices, as.integer),
  target = lapply(seq_len(n_spots), \(i) {
    sort(as.integer(which(hex_dist[i, ] <= 110) - 1L))
  }),
  info = "sp graph - the radius layout keeps every pair within the cut-off"
)

expect_equal(
  current = lapply(radius_graph$indices, as.integer),
  target = lapply(hex_binary$indices, as.integer),
  info = "sp graph - a 110 px radius reproduces the honeycomb exactly"
)

## weightings ------------------------------------------------------------------

expect_true(
  current = all(unlist(hex_binary$weights) == 1),
  info = "sp graph - binary weighting gives weights of one"
)

inv_dist <- build(
  hex_coords,
  params_sp_graph("radius", "inverse_distance", radius = 110, power = 2)
)

expect_equal(
  current = as.numeric(unlist(inv_dist$weights)),
  target = as.numeric(unlist(lapply(seq_len(n_spots), \(i) {
    hex_dist[i, inv_dist$indices[[i]] + 1L]^-2
  }))),
  tolerance = 1e-6,
  info = "sp graph - inverse distance weighting is d^-power"
)

gaussian <- build(
  hex_coords,
  params_sp_graph("radius", "gaussian", radius = 110, bandwidth = 50)
)

expect_equal(
  current = as.numeric(unlist(gaussian$weights)),
  target = as.numeric(unlist(lapply(seq_len(n_spots), \(i) {
    exp(-hex_dist[i, gaussian$indices[[i]] + 1L]^2 / (2 * 50^2))
  }))),
  tolerance = 1e-6,
  info = "sp graph - gaussian weighting is exp(-d^2 / (2 b^2))"
)

# a NULL bandwidth resolves to the median nearest-neighbour distance inside
# Rust, so the nearest edge of every spot lands on exp(-1/2)
gaussian_auto <- build(
  hex_coords,
  params_sp_graph("radius", "gaussian", radius = 110)
)
nearest_weight <- vapply(
  seq_len(n_spots),
  \(i) max(gaussian_auto$weights[[i]]),
  numeric(1)
)

expect_equal(
  current = max(nearest_weight),
  target = exp(-0.5),
  tolerance = 1e-5,
  info = "sp graph - a NULL bandwidth resolves to the median NN distance"
)

row_norm <- build(
  hex_coords,
  params_sp_graph("hex", "binary", row_normalise = TRUE),
  hex_grid
)

expect_equal(
  current = vapply(row_norm$weights, sum, numeric(1)),
  target = rep(1, n_spots),
  tolerance = 1e-6,
  info = "sp graph - row normalisation makes every row sum to one"
)

expect_equal(
  current = as.numeric(stats::median(row_norm$degree)),
  target = 2,
  tolerance = 1e-6,
  info = "sp graph - row normalised degree is 2, not the neighbour count"
)

## failure modes ---------------------------------------------------------------

expect_error(
  current = rs_sp_build_graph(
    coords = hex_coords,
    graph_params = params_sp_graph("hex", "binary"),
    knn_params = NULL,
    seed = 42L,
    verbose = 0L
  ),
  info = "sp graph - a lattice layout without array coordinates errors"
)

expect_error(
  current = rs_sp_build_graph(
    coords = hex_coords,
    graph_params = with_array(
      params_sp_graph("gaussian", "gaussian"),
      hex_grid
    ),
    knn_params = NULL,
    seed = 42L,
    verbose = 0L
  ),
  info = "sp graph - an unknown layout errors"
)

## through the class -----------------------------------------------------------

fixture <- sp_make_visium(
  file.path(test_temp_dir, "sample_a"),
  layout = "hex",
  n_rows = 10L,
  n_cols = 10L,
  n_genes = 30L
)
object <- sp_load_visium_fixture(fixture, file.path(test_temp_dir, "out_a"))

obs <- sp_obs_in_graph_order(object, "sample_a")
coords <- get_spatial_coords(object, "sample_a")

# The sp_cache getters and setters are S3 generics pointed at S7 classes, so
# they only fire when registered through S7::method(). A plain
# `get_per_sample_spatial_graph.SpatialSpot` would silently do nothing and this
# would come back non-NULL forever.
expect_null(
  current = get_per_sample_spatial_graph(object, "sample_a"),
  info = "sp graph - nothing is cached before the graph is built"
)

object <- build_spatial_graph_sp(object, .verbose = FALSE)
cached <- get_per_sample_spatial_graph(object, "sample_a")

expect_true(
  current = inherits(cached, "CsparseMatrix"),
  info = "sp graph - the built graph is cached on a SpatialSpot"
)

direct <- rs_sp_build_graph(
  coords = coords,
  graph_params = with_array(params_sp_graph(), obs),
  knn_params = NULL,
  seed = 42L,
  verbose = 0L
)

expect_equal(
  current = bixverse:::.sp_sparse_to_adjacency(cached),
  target = list(
    indices = lapply(direct$indices, as.integer),
    weights = lapply(direct$weights, as.numeric)
  ),
  info = "sp graph - the cached matrix round trips to the direct Rust call"
)

expect_equal(
  current = lapply(direct$indices, as.integer),
  target = sp_lattice_neighbours(obs$array_row, obs$array_col, "hex"),
  info = "sp graph - the method injects the array coordinates in obs order"
)

# the coordinates carry no lattice information, so a graph that ignored the
# injected array coordinates would come out different
knn_object <- build_spatial_graph_sp(
  object,
  graph_params = params_sp_graph("knn", k = 4L),
  .verbose = FALSE
)

expect_equal(
  current = lengths(
    bixverse:::.sp_sparse_to_adjacency(
      get_per_sample_spatial_graph(knn_object, "sample_a")
    )$indices
  ),
  target = rep(4L, fixture$n_spots),
  info = "sp graph - the layout choice reaches Rust through the method"
)

expect_error(
  current = build_spatial_graph_sp(
    object,
    exp_id = "not_a_sample",
    .verbose = FALSE
  ),
  info = "sp graph - an unknown exp_id errors"
)

## a subset agrees with its parent ---------------------------------------------

subset_graph <- get_per_sample_spatial_graph(
  build_spatial_graph_sp(
    SpatialSpotSubset(object, exp_id = "sample_a"),
    .verbose = FALSE
  ),
  "sample_a"
)

expect_equal(
  current = as.matrix(subset_graph),
  target = as.matrix(cached),
  info = "sp graph - a subset builds the same graph as its parent"
)

## several samples -------------------------------------------------------------

# The graph is per sample, so two samples of different shape have to end up
# with two adjacencies of different size rather than one shared one.

fixture_b <- sp_make_visium(
  file.path(test_temp_dir, "sample_b"),
  layout = "hex",
  n_rows = 6L,
  n_cols = 6L,
  n_genes = 30L,
  seed = 77L
)

multi_dir <- file.path(test_temp_dir, "out_multi")
unlink(multi_dir, recursive = TRUE)
dir.create(multi_dir, recursive = TRUE, showWarnings = FALSE)

prescan <- prescan_visium_dirs(
  dirs = c(fixture$visium_dir, fixture_b$visium_dir),
  exp_ids = c("sample_a", "sample_b"),
  .verbose = FALSE
)

multi <- load_multi_visium(
  SpatialSpot(dir_data = multi_dir),
  prescan_result = prescan,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 0L,
    min_lib_size = 0L,
    min_cells = 0L
  ),
  .verbose = FALSE
)

multi <- build_spatial_graph_sp(multi, .verbose = FALSE)

expect_equal(
  current = dim(get_per_sample_spatial_graph(multi, "sample_a")),
  target = rep(fixture$n_spots, 2L),
  info = "sp graph - the first sample gets a graph of its own size"
)

expect_equal(
  current = dim(get_per_sample_spatial_graph(multi, "sample_b")),
  target = rep(fixture_b$n_spots, 2L),
  info = "sp graph - the second one gets a smaller, separate graph"
)

expect_equal(
  current = as.matrix(get_per_sample_spatial_graph(multi, "sample_a")),
  target = as.matrix(cached),
  info = "sp graph - a sample's graph does not change by loading a second one"
)

# restricting to one exp_id must leave the other sample untouched
partial <- build_spatial_graph_sp(
  multi,
  exp_id = "sample_b",
  graph_params = params_sp_graph("knn", k = 3L),
  .verbose = FALSE
)

expect_equal(
  current = unique(lengths(
    bixverse:::.sp_sparse_to_adjacency(
      get_per_sample_spatial_graph(partial, "sample_b")
    )$indices
  )),
  target = 3L,
  info = "sp graph - exp_id restricts the rebuild to that sample"
)

expect_equal(
  current = as.matrix(get_per_sample_spatial_graph(partial, "sample_a")),
  target = as.matrix(cached),
  info = "sp graph - and leaves the other sample's graph alone"
)
