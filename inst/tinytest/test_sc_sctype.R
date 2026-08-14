# sctype tests -----------------------------------------------------------------

source("helper_sc.R", local = TRUE)

set.seed(123L)

test_temp_dir <- sc_test_dir("sctype")

## testing parameters ----------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
# hvg / pca
hvg_to_keep <- 30L
no_pcs <- 15L

## synthetic test data ---------------------------------------------------------

# the default configuration makes genes 1:10 markers of cell_type_1, 11:20 of
# cell_type_2 and 21:30 of cell_type_3
single_cell_test_data <- sc_test_fixture(
  min_lib_size = min_lib_size,
  min_genes_exp = min_genes_exp,
  min_cells_exp = min_cells_exp,
  hvg_to_keep = hvg_to_keep,
  no_pcs = no_pcs
)

marker_df <- data.table::data.table(
  cell_type = rep(sprintf("cell_type_%i", 1:3), each = 10),
  gene_id = sprintf("gene_%03d", 1:30)
)

sc_qc_param <- sc_test_qc_params(single_cell_test_data, target_size = 1000)

load_sc_object <- function(name) {
  sc_test_object(
    sc_test_dir(name, test_temp_dir),
    single_cell_test_data,
    sc_qc_param = sc_qc_param
  )
}

## processed object ------------------------------------------------------------

sc_object <- load_sc_object("processed")

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = no_pcs,
  pca_params = params_sc_pca(),
  .verbose = FALSE
)

sc_object <- find_neighbours_sc(sc_object, .verbose = FALSE)
sc_object <- find_clusters_sc(sc_object)

n_cells <- length(get_cells_to_keep(sc_object))
cell_grps <- unlist(sc_object[["cell_grp"]])
cluster_labels <- unlist(sc_object[["leiden_clustering"]])

# tests ------------------------------------------------------------------------

## marker preparation ----------------------------------------------------------

cell_markers <- prepare_cell_markers(
  obj = sc_object,
  marker_df = data.table::copy(marker_df)
)

expect_equal(
  current = names(cell_markers),
  target = sprintf("cell_type_%i", 1:3),
  info = "prepare_cell_markers - one entry per cell type"
)

## scoring ---------------------------------------------------------------------

sc_type_res <- calc_sc_type_scores(
  object = sc_object,
  cell_marker_list = cell_markers,
  .verbose = FALSE
)

expect_true(
  current = inherits(sc_type_res, "ScTypeResults"),
  info = "calc_sc_type_scores - returns an ScTypeResults"
)

expect_equal(
  current = sc_type_res$n_cells,
  target = n_cells,
  info = "calc_sc_type_scores - scores every kept cell"
)

expect_equal(
  current = dim(get_scores(sc_type_res)),
  target = c(n_cells, 3L),
  info = "calc_sc_type_scores - score matrix is cells x cell types"
)

## per-cell assignment ---------------------------------------------------------

cell_res <- assign_sc_type(
  object = sc_object,
  sc_type_res = sc_type_res,
  .verbose = FALSE
)

expect_true(
  current = inherits(cell_res, "ScTypeCellResults"),
  info = "assign_sc_type - returns an ScTypeCellResults"
)

expect_equal(
  current = length(cell_res$assignments),
  target = n_cells,
  info = "assign_sc_type - one call per cell"
)

expect_true(
  current = all(
    stats::na.omit(cell_res$assignments) %in% cell_res$cell_types
  ),
  info = "assign_sc_type - every call is a scored cell type"
)

expect_true(
  current = !is.null(cell_res$agreement) &&
    checkmate::qtest(cell_res$agreement, sprintf("N%i[0,1]", n_cells)),
  info = "assign_sc_type - neighbour agreement returned with an sNN graph"
)

expect_null(
  current = cell_res$composition,
  info = "assign_sc_type - no composition without a cluster column"
)

expect_true(
  current = mean(cell_res$assignments == cell_grps, na.rm = TRUE) > 0.9,
  info = "assign_sc_type - per-cell calls recover the cell groups"
)

### reproducibility ------------------------------------------------------------

cell_res_rerun <- assign_sc_type(
  object = sc_object,
  sc_type_res = sc_type_res,
  .verbose = FALSE
)

expect_equal(
  current = cell_res_rerun$assignments,
  target = cell_res$assignments,
  info = "assign_sc_type - identical calls on a rerun"
)

### no smoothing ---------------------------------------------------------------

sc_object_raw <- load_sc_object("no_graph")

expect_warning(
  current = assign_sc_type(
    object = sc_object_raw,
    sc_type_res = sc_type_res,
    .verbose = FALSE
  ),
  info = "assign_sc_type - warns when no sNN graph is present"
)

cell_res_raw <- suppressWarnings(assign_sc_type(
  object = sc_object_raw,
  sc_type_res = sc_type_res,
  .verbose = FALSE
))

expect_null(
  current = cell_res_raw$agreement,
  info = "assign_sc_type - no agreement without a graph"
)

## composition and hybrid ------------------------------------------------------

hybrid_res <- assign_sc_type(
  object = sc_object,
  sc_type_res = sc_type_res,
  cluster_col = "leiden_clustering",
  .verbose = FALSE
)

comp <- hybrid_res$composition

expect_true(
  current = data.table::is.data.table(comp),
  info = "assign_sc_type - composition is a data.table"
)

expect_equal(
  current = nrow(comp),
  target = length(unique(cluster_labels)),
  info = "assign_sc_type - one composition row per non-empty cluster"
)

expect_true(
  current = checkmate::qtest(comp$purity, sprintf("N%i[0,1]", nrow(comp))),
  info = "assign_sc_type - purity is a fraction"
)

expect_equal(
  current = unname(rowSums(hybrid_res$counts)) + comp$n_unknown,
  target = comp$n_cells,
  info = "assign_sc_type - counts plus unknowns add up per cluster"
)

expect_equal(
  current = comp$n_cells,
  target = as.integer(table(cluster_labels)[as.character(comp$cluster_id)]),
  info = "assign_sc_type - cluster sizes match the cluster labels"
)

expect_true(
  current = all(
    tapply(
      hybrid_res$hybrid[
        !comp$cluster_mixed[
          match(cluster_labels, comp$cluster_id)
        ]
      ],
      cluster_labels[
        !comp$cluster_mixed[
          match(cluster_labels, comp$cluster_id)
        ]
      ],
      \(x) length(unique(x)) == 1
    )
  ),
  info = "assign_sc_type - pure clusters get one blanket label"
)

### threshold behaviour --------------------------------------------------------

hybrid_blanket <- assign_sc_type(
  object = sc_object,
  sc_type_res = sc_type_res,
  cluster_col = "leiden_clustering",
  sctype_cell_params = params_sctype_cells(purity_threshold = 0),
  .verbose = FALSE
)

expect_true(
  current = !any(hybrid_blanket$composition$cluster_mixed),
  info = "assign_sc_type - purity_threshold of 0 keeps every cluster call"
)

expect_true(
  current = all(tapply(
    hybrid_blanket$hybrid,
    cluster_labels,
    \(x) length(unique(x)) == 1
  )),
  info = "assign_sc_type - purity_threshold of 0 gives uniform cluster labels"
)

hybrid_split <- assign_sc_type(
  object = sc_object,
  sc_type_res = sc_type_res,
  cluster_col = "leiden_clustering",
  sctype_cell_params = params_sctype_cells(purity_threshold = 1),
  .verbose = FALSE
)

expect_equal(
  current = hybrid_split$hybrid,
  target = hybrid_split$assignments,
  info = "assign_sc_type - purity_threshold of 1 falls back to per-cell calls"
)

## calibration -----------------------------------------------------------------

cell_res_z <- assign_sc_type(
  object = sc_object,
  sc_type_res = sc_type_res,
  sctype_cell_params = params_sctype_cells(calibration = "column_z"),
  .verbose = FALSE
)

expect_equal(
  current = length(cell_res_z$assignments),
  target = n_cells,
  info = "assign_sc_type - column_z calibration runs"
)

## parameter validation --------------------------------------------------------

expect_error(
  current = params_sctype_cells(alpha = 2),
  info = "params_sctype_cells - alpha needs to be in [0, 1]"
)

expect_error(
  current = params_sctype_cells(calibration = "nope"),
  info = "params_sctype_cells - unknown calibration is rejected"
)

expect_error(
  current = params_sctype_cells(tolerance = 0),
  info = "params_sctype_cells - tolerance needs to be positive"
)

expect_error(
  current = assign_sc_type(
    object = sc_object,
    sc_type_res = sc_type_res,
    sctype_cell_params = list(alpha = 0.5),
    .verbose = FALSE
  ),
  info = "assign_sc_type - incomplete parameter list is rejected"
)

# clean up ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir)
