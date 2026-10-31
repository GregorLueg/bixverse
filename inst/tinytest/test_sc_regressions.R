# regressions fixed in 0.9.0 ---------------------------------------------------

sc <- demo_single_cells()
dir_data <- S7::prop(sc, "dir_data")

## group codes -----------------------------------------------------------------

expect_equal(
  current = bixverse:::.group_codes(c("b", "B", "a", "b"), "x"),
  target = c(2L, 0L, 1L, 2L),
  info = "group codes - code point order, B before a before b"
)

expect_error(
  current = bixverse:::.group_codes(c("a", NA), "batch"),
  pattern = "`batch` has 1 NA labels",
  info = "group codes - NA labels error"
)

## params ----------------------------------------------------------------------

milo <- params_sc_miloR(knn = list(m = 8L, ef_search = 50L))

expect_equal(
  current = c(milo$m, milo$ef_search),
  target = c(8L, 50L),
  info = "params_sc_miloR - kNN overrides reach the list"
)

expect_false(
  current = any(c("nn_max_iter", "rho") %in% names(milo)),
  info = "params_sc_miloR - no dead kNN keys"
)

expect_error(
  current = params_sc_miloR(index_type = "ivf"),
  info = "params_sc_miloR - ivf has no Milo index"
)

expect_error(
  current = params_sc_miloR(knn = list(not_a_key = 1L)),
  info = "params_sc_miloR - unknown kNN keys error"
)

expect_false(
  current = "nn_max_iter" %in% names(params_sc_vision()$knn),
  info = "params_sc_vision - no dead nn_max_iter"
)

expect_error(
  current = params_sc_harmony(sigma = 0),
  info = "params_sc_harmony - sigma must be positive"
)

expect_false(
  current = "randomised" %in% names(params_sc_pca()),
  info = "params_sc_pca - deprecated randomised is gone"
)

expect_error(
  current = generate_seacells_sc(sc, .verbose = FALSE),
  info = "generate_seacells_sc - seacell_params is required"
)

## kNN -------------------------------------------------------------------------

knn_3 <- generate_knn_sc(
  sc,
  no_embd_to_use = 3L,
  neighbours_params = params_sc_neighbours(
    knn = list(knn_method = "exhaustive")
  ),
  .verbose = FALSE
)
ref_3 <- rs_sc_knn(
  embd = get_pca_factors(sc)[, 1:3],
  knn_params = list(knn_method = "exhaustive", ann_dist = "euclidean"),
  validate_index = FALSE,
  verbose = 0L,
  seed = 42L
)

expect_equal(
  current = get_knn_mat(knn_3),
  target = ref_3,
  info = "generate_knn_sc - honours no_embd_to_use"
)

# row 0: neighbours stored by column index with descending distances
bbknn <- rs_bbknn_filtering(
  indptr = c(0L, 3L, 5L),
  indices = c(0L, 1L, 2L, 0L, 1L),
  data = c(3, 2, 1, 0.5, 0.5),
  no_neighbours_to_keep = 2L
)

expect_equal(
  current = bbknn$indices[1, ],
  target = c(2L, 1L),
  info = "rs_bbknn_filtering - keeps the nearest, not the lowest index"
)

## clustering ------------------------------------------------------------------

sc <- find_clusters_sc(sc, cluster_algorithm = "louvain", .verbose = FALSE)

expect_true(
  current = "louvain_clustering" %in% names(get_sc_obs(sc)),
  info = "find_clusters_sc - takes .verbose, louvain gets its own column"
)

## quality control -------------------------------------------------------------

expect_error(
  current = run_cell_qc(
    metrics = list(lib_size = c(1, 2, 3, 100)),
    cells_to_keep = 0:2
  ),
  pattern = "cells_to_keep",
  info = "run_cell_qc - cells_to_keep length is checked"
)

expect_error(
  current = run_cell_qc(
    metrics = list(lib_size = c(1, 2, 3, 100)),
    cells_to_keep = 0:3,
    directions = c(lib_size = "sideways")
  ),
  info = "run_cell_qc - unknown direction errors"
)

## markers ---------------------------------------------------------------------

markers <- data.table::data.table(
  gene_id = get_gene_names(sc)[1:4],
  cell_type = c("a", "a", "b", "b")
)
invisible(prepare_cell_markers(obj = sc, marker_df = markers))

expect_false(
  current = "gene_idx" %in% names(markers),
  info = "prepare_cell_markers - leaves the caller's table alone"
)

## pipelines -------------------------------------------------------------------

expect_error(
  current = bixverse:::new_sc_step(
    "pca",
    calculate_pca_sc,
    list(no_pc = 5L)
  ),
  pattern = "unknown argument",
  info = "pipeline - a misspelt step argument fails at build time"
)

sc <- set_sc_new_obs_col(
  sc,
  col_name = "num_grp",
  new_data = as.integer(factor(get_sc_obs(sc, filtered = TRUE)$cell_grp))
)

per_batch <- apply_pipeline_per_group(
  sc_pipeline() %>>% step_hvg_sc(hvg_no = 20L, .verbose = FALSE),
  sc,
  group_col = "num_grp"
)

expect_true(
  current = length(per_batch) > 1L,
  info = "pipeline - per group works on a numeric column"
)

unlink(dir_data, recursive = TRUE, force = TRUE)
