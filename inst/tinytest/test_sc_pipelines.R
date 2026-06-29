# sc pipeline ------------------------------------------------------------------

set.seed(123L)

test_temp_dir <- file.path(tempdir(), "sc_pipeline")
dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

## parameters ------------------------------------------------------------------

min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
hvg_to_keep <- 30L
no_pcs <- 15L

## set up parent object --------------------------------------------------------

single_cell_test_data <- generate_single_cell_test_data()

sc_qc_param <- params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp,
  target_size = 1000
)

sc_object <- SingleCells(dir_data = test_temp_dir)

sc_object <- load_r_data(
  object = sc_object,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  var = single_cell_test_data$var,
  sc_qc_param = sc_qc_param,
  streaming = 0L,
  .verbose = FALSE
)

# tests ------------------------------------------------------------------------

## construction ----------------------------------------------------------------

empty <- sc_pipeline()

expect_true(
  current = inherits(empty, "ScPipeline"),
  info = "sc_pipeline() returns a ScPipeline"
)

expect_equal(
  current = length(empty),
  target = 0L,
  info = "empty pipeline has length 0"
)

s_hvg <- step_hvg_sc(hvg_no = hvg_to_keep, .verbose = FALSE)

expect_true(
  current = inherits(s_hvg, "ScStep"),
  info = "step_hvg_sc() returns a ScStep"
)

expect_equal(
  current = s_hvg$name,
  target = "hvg",
  info = "step carries human-readable name"
)

expect_true(
  current = is.function(s_hvg$fn),
  info = "step carries the underlying function"
)

expect_equal(
  current = s_hvg$args$hvg_no,
  target = hvg_to_keep,
  info = "step captures named args verbatim"
)

## chaining --------------------------------------------------------------------

### pipeline %>>% step ---------------------------------------------------------

p1 <- sc_pipeline() %>>% step_hvg_sc(hvg_no = hvg_to_keep, .verbose = FALSE)

expect_true(
  current = inherits(p1, "ScPipeline"),
  info = "pipeline %>>% step returns ScPipeline"
)

expect_equal(
  current = length(p1),
  target = 1L,
  info = "appending one step gives length 1"
)

### step %>>% step -------------------------------------------------------------

p2 <- step_hvg_sc(hvg_no = hvg_to_keep, .verbose = FALSE) %>>%
  step_pca_sc(no_pcs = no_pcs, .verbose = FALSE)

expect_true(
  current = inherits(p2, "ScPipeline"),
  info = "step %>>% step wraps into a ScPipeline"
)

expect_equal(
  current = length(p2),
  target = 2L,
  info = "step %>>% step has 2 steps"
)

expect_equal(
  current = vapply(p2$steps, `[[`, character(1), "name"),
  target = c("hvg", "pca"),
  info = "step order preserved through %>>%"
)

### longer chain --------------------------------------------------------------

p3 <- sc_pipeline() %>>%
  step_hvg_sc(hvg_no = hvg_to_keep, .verbose = FALSE) %>>%
  step_pca_sc(no_pcs = no_pcs, .verbose = FALSE) %>>%
  step_neighbours_sc(.verbose = FALSE) %>>%
  step_clusters_sc(name = "leiden_clustering")

expect_equal(
  current = length(p3),
  target = 4L,
  info = "4-step chain has length 4"
)

expect_equal(
  current = vapply(p3$steps, `[[`, character(1), "name"),
  target = c("hvg", "pca", "neighbours", "clusters"),
  info = "4-step chain preserves order"
)

### rhs must be a step --------------------------------------------------------

expect_error(
  current = sc_pipeline() %>>% "not_a_step",
  info = "%>>% errors when rhs is not a ScStep"
)

## primitives ------------------------------------------------------------------

expect_silent(
  current = capture.output(print(empty)),
  info = "print() works on empty pipeline"
)

expect_silent(
  current = capture.output(print(p3)),
  info = "print() works on populated pipeline"
)

expect_silent(
  current = capture.output(print(s_hvg)),
  info = "print() works on a ScStep"
)

## apply_pipeline --------------------------------------------------------------

### empty pipeline returns object unchanged -----------------------------------

obj_before <- sc_object
obj_after <- apply_pipeline(sc_pipeline(), sc_object)

expect_equal(
  current = get_hvg(obj_after),
  target = get_hvg(obj_before),
  info = "empty pipeline leaves object unchanged"
)

### single step runs underlying function --------------------------------------

p_hvg <- sc_pipeline() %>>%
  step_hvg_sc(hvg_no = hvg_to_keep, .verbose = FALSE)

sc_object <- apply_pipeline(p_hvg, sc_object)

expect_true(
  current = length(get_hvg(sc_object)) == hvg_to_keep,
  info = "single-step pipeline ran HVG and set indices"
)

### multi-step runs each step in order ----------------------------------------

p_hvg_pca <- sc_pipeline() %>>%
  step_hvg_sc(hvg_no = hvg_to_keep, .verbose = FALSE) %>>%
  step_pca_sc(no_pcs = no_pcs, .verbose = FALSE)

sc_object_2 <- apply_pipeline(p_hvg_pca, SingleCells(dir_data = test_temp_dir))
# reload counts on the fresh handle
sc_object_2 <- load_r_data(
  object = sc_object_2,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  var = single_cell_test_data$var,
  sc_qc_param = sc_qc_param,
  streaming = 0L,
  .verbose = FALSE
)
sc_object_2 <- apply_pipeline(p_hvg_pca, sc_object_2)

expect_true(
  current = length(get_hvg(sc_object_2)) == hvg_to_keep,
  info = "multi-step: HVG set"
)

expect_equal(
  current = ncol(get_pca_factors(sc_object_2)),
  target = no_pcs,
  info = "multi-step: PCA factors set with correct dimensionality"
)

### wrong order surfaces the underlying warning -------------------------------

# PCA before HVG -> calculate_pca_sc warns about missing HVGs
p_wrong <- sc_pipeline() %>>%
  step_pca_sc(no_pcs = no_pcs, .verbose = FALSE) %>>%
  step_hvg_sc(hvg_no = hvg_to_keep, .verbose = FALSE)

fresh <- SingleCells(dir_data = file.path(tempdir(), "sc_pipeline_wrong"))
dir.create(dirname(fresh@dir_data), recursive = TRUE, showWarnings = FALSE)
fresh <- load_r_data(
  object = fresh,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  var = single_cell_test_data$var,
  sc_qc_param = sc_qc_param,
  streaming = 0L,
  .verbose = FALSE
)

expect_warning(
  current = apply_pipeline(p_wrong, fresh),
  info = "underlying warnings propagate through apply_pipeline"
)

### works on a subset ---------------------------------------------------------

subset_obj <- SingleCellsSubset(
  sc_object = sc_object,
  grouping_column = "cell_grp",
  group = "cell_type_1"
)

p_sub <- sc_pipeline() %>>%
  step_hvg_sc(hvg_no = hvg_to_keep, .verbose = FALSE) %>>%
  step_pca_sc(no_pcs = no_pcs, .verbose = FALSE)

subset_obj <- apply_pipeline(p_sub, subset_obj)

expect_true(
  current = length(get_hvg(subset_obj)) == hvg_to_keep,
  info = "pipeline dispatches on SingleCellsSubset (HVG)"
)

expect_equal(
  current = ncol(get_pca_factors(subset_obj)),
  target = no_pcs,
  info = "pipeline dispatches on SingleCellsSubset (PCA)"
)

## apply_pipeline_per_group ----------------------------------------------------

p_group <- sc_pipeline() %>>%
  step_hvg_sc(hvg_no = hvg_to_keep, .verbose = FALSE) %>>%
  step_pca_sc(no_pcs = no_pcs, .verbose = FALSE)

per_group <- apply_pipeline_per_group(
  pipeline = p_group,
  object = sc_object,
  group_col = "cell_grp"
)

expected_groups <- unique(as.character(
  get_sc_obs(sc_object, filtered = TRUE)$cell_grp
))

expect_true(
  current = setequal(names(per_group), expected_groups),
  info = "per-group apply covers all groups"
)

expect_true(
  current = all(vapply(
    per_group,
    function(x) S7::S7_inherits(x, SingleCellsSubset),
    logical(1)
  )),
  info = "per-group apply returns SingleCellsSubset objects"
)

expect_true(
  current = all(vapply(
    per_group,
    function(x) length(get_hvg(x)) == hvg_to_keep,
    logical(1)
  )),
  info = "pipeline ran on every group (HVG set in each subset)"
)

### restricting groups --------------------------------------------------------

restricted <- apply_pipeline_per_group(
  pipeline = p_group,
  object = sc_object,
  group_col = "cell_grp",
  groups = expected_groups[1]
)

expect_equal(
  current = names(restricted),
  target = expected_groups[1],
  info = "groups argument restricts output"
)

### bad inputs ----------------------------------------------------------------

expect_error(
  current = apply_pipeline_per_group(
    pipeline = p_group,
    object = sc_object,
    group_col = "not_a_column"
  ),
  info = "errors on missing grouping column"
)

## teardown --------------------------------------------------------------------

unlink(test_temp_dir, recursive = TRUE, force = TRUE)
unlink(
  file.path(tempdir(), "sc_pipeline_wrong"),
  recursive = TRUE,
  force = TRUE
)
