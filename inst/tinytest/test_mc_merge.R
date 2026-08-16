# meta cell merging ------------------------------------------------------------

source("helper_sc.R", local = TRUE)

set.seed(123L)

test_temp_dir <- sc_test_dir("mc_merge")

## namespace regression --------------------------------------------------------

# `Matrix` is only ever called as `Matrix::`. Without a NAMESPACE import its
# namespace stays unloaded, and the first S4 dispatch on a `dgRMatrix` that came
# out of a readRDS/qs2 load misfires: `nrow()` gives NULL, `[` refuses to
# subset, and Rust reports "nrow missing or not a non-negative whole number".
# Cannot be reproduced in-process, since Matrix is loaded by the time this runs.
expect_true(
  current = "Matrix" %in% names(getNamespaceImports("bixverse")),
  info = "mc merge - bixverse imports from Matrix so its namespace loads eagerly"
)

## hand-built inputs -----------------------------------------------------------

# Small `MetaCells` objects with a gene space we control, so the remapping
# branches of `.merge_mc_counts()` can be driven directly. `dense` is the same
# matrix as a base one, used to build the reference for the merged counts.
build_mc <- function(n_mc, gene_ids, seed, n_cells = 100L, id_prefix = NULL) {
  set.seed(seed)
  m <- Matrix::rsparsematrix(
    n_mc,
    length(gene_ids),
    density = 0.5,
    repr = "R"
  )
  m@x <- abs(round(m@x * 10)) + 1

  object <- MetaCells(
    meta_cell_data = list(
      aggregated = list(
        indptr = m@p,
        indices = m@j,
        raw_counts = m@x,
        norm_counts = m@x / 2,
        nrow = n_mc,
        ncol = length(gene_ids)
      ),
      assignments = list(
        assignments = as.integer(rep(seq_len(n_mc), length.out = n_cells)),
        metacells = purrr::map(seq_len(n_mc), \(i) as.integer(seq_len(3L) + i)),
        n_cells = n_cells,
        n_metacells = n_mc,
        n_unassigned = 0L
      )
    ),
    var_data = data.table::data.table(
      gene_idx = seq_along(gene_ids),
      gene_id = gene_ids
    ),
    meta_cell_method = "meta_cells_hdwgcna",
    obs_ids = if (is.null(id_prefix)) {
      NULL
    } else {
      sprintf("%s_%02d", id_prefix, seq_len(n_mc))
    }
  )

  list(object = object, dense = as.matrix(m), genes = gene_ids)
}

# Row-bind the dense counterparts onto `target` the slow, obvious way.
dense_reference <- function(parts, target) {
  do.call(
    rbind,
    purrr::map(parts, \(p) {
      out <- matrix(0, nrow = nrow(p$dense), ncol = length(target))
      hit <- match(target, p$genes)
      out[, !is.na(hit)] <- p$dense[, hit[!is.na(hit)], drop = FALSE]
      out
    })
  )
}

genes_a <- paste0("gene_", 1:12)
# reversed order, one gene of `a` missing, one gene of its own
genes_b <- paste0("gene_", c(11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 13))

part_a <- build_mc(6L, genes_a, seed = 1L)
part_b <- build_mc(5L, genes_b, seed = 2L)

# tests ------------------------------------------------------------------------

## gene space ------------------------------------------------------------------

### intersect ------------------------------------------------------------------

merged_int <- suppressWarnings(merge_meta_cells(
  inputs = list(a = part_a$object, b = part_b$object),
  feature_space = "intersect",
  .verbose = FALSE
))

expect_equal(
  current = S7::prop(merged_int, "var_table")$gene_id,
  target = intersect(genes_a, genes_b),
  info = "mc merge - intersect keeps the shared genes in first-input order"
)

expect_equal(
  current = S7::prop(merged_int, "dims"),
  target = c(11L, 11L),
  info = "mc merge - intersect dims are (sum of meta cells, shared genes)"
)

expect_equal(
  current = unname(as.matrix(S7::prop(merged_int, "data")$raw)),
  target = dense_reference(
    list(part_a, part_b),
    S7::prop(merged_int, "var_table")$gene_id
  ),
  info = "mc merge - intersect raw counts match a dense row-bind"
)

expect_equal(
  current = unname(as.matrix(S7::prop(merged_int, "data")$norm)),
  target = dense_reference(
    list(part_a, part_b),
    S7::prop(merged_int, "var_table")$gene_id
  ) /
    2,
  info = "mc merge - intersect norm counts ride the same index structure"
)

expect_true(
  current = validObject(S7::prop(merged_int, "data")$raw),
  info = "mc merge - intersect produces a valid sparse matrix"
)

### union ----------------------------------------------------------------------

merged_uni <- suppressWarnings(merge_meta_cells(
  inputs = list(a = part_a$object, b = part_b$object),
  feature_space = "union",
  .verbose = FALSE
))

expect_equal(
  current = S7::prop(merged_uni, "var_table")$gene_id,
  target = c(genes_a, setdiff(genes_b, genes_a)),
  info = "mc merge - union appends the genes the first input does not hold"
)

expect_equal(
  current = unname(as.matrix(S7::prop(merged_uni, "data")$raw)),
  target = dense_reference(
    list(part_a, part_b),
    S7::prop(merged_uni, "var_table")$gene_id
  ),
  info = "mc merge - union raw counts match a dense row-bind"
)

expect_true(
  current = validObject(S7::prop(merged_uni, "data")$raw),
  info = "mc merge - union produces a valid sparse matrix"
)

# `gene_13` only exists in `b`, `gene_12` only in `a`
union_raw <- as.matrix(S7::prop(merged_uni, "data")$raw)
union_genes <- S7::prop(merged_uni, "var_table")$gene_id

expect_true(
  current = all(union_raw[1:6, which(union_genes == "gene_13")] == 0),
  info = "mc merge - union leaves structural zeros for genes absent from a"
)

expect_true(
  current = all(union_raw[7:11, which(union_genes == "gene_12")] == 0),
  info = "mc merge - union leaves structural zeros for genes absent from b"
)

### degenerate cases -----------------------------------------------------------

merged_one <- merge_meta_cells(
  inputs = list(only = part_a$object),
  .verbose = FALSE
)

expect_equal(
  current = unname(as.matrix(S7::prop(merged_one, "data")$raw)),
  target = part_a$dense,
  info = "mc merge - a single input round trips unchanged"
)

expect_error(
  current = merge_meta_cells(
    inputs = list(
      a = part_a$object,
      z = build_mc(3L, paste0("other_", 1:4), seed = 3L)$object
    ),
    .verbose = FALSE
  ),
  info = "mc merge - an empty gene intersection errors"
)

## observation table -----------------------------------------------------------

obs_merged <- S7::prop(merged_int, "obs_table")

expect_true(
  current = "source_id" %in% names(obs_merged),
  info = "mc merge - the observation table gains a source_id column"
)

expect_equal(
  current = names(obs_merged)[1:3],
  target = c("meta_cell_idx", "meta_cell_id", "source_id"),
  info = "mc merge - the identifying columns come first"
)

expect_equal(
  current = obs_merged$meta_cell_idx,
  target = 1:11,
  info = "mc merge - meta_cell_idx is renumbered over the merged object"
)

expect_equal(
  current = obs_merged$source_id,
  target = c(rep("a", 6), rep("b", 5)),
  info = "mc merge - source_id follows the input order"
)

expect_true(
  current = all(grepl("^(a|b)__meta_cell_", obs_merged$meta_cell_id)),
  info = "mc merge - meta cell identifiers are prefixed with their source"
)

expect_equal(
  current = rownames(S7::prop(merged_int, "data")$raw),
  target = obs_merged$meta_cell_id,
  info = "mc merge - the count matrix carries the merged identifiers"
)

expect_equal(
  current = obs_merged$original_cell_idx,
  target = c(
    S7::prop(part_a$object, "obs_table")$original_cell_idx,
    S7::prop(part_b$object, "obs_table")$original_cell_idx
  ),
  info = "mc merge - original_cell_idx is carried over untouched"
)

### identifiers ----------------------------------------------------------------

expect_error(
  current = suppressWarnings(merge_meta_cells(
    inputs = list(a = part_a$object, b = part_b$object),
    prefix_ids = FALSE,
    .verbose = FALSE
  )),
  info = "mc merge - duplicated identifiers error when not prefixing"
)

merged_unprefixed <- merge_meta_cells(
  inputs = list(
    a = build_mc(6L, genes_a, seed = 1L, id_prefix = "left")$object,
    b = build_mc(5L, genes_a, seed = 4L, id_prefix = "right")$object
  ),
  source_ids = c("first", "second"),
  prefix_ids = FALSE,
  .verbose = FALSE
)

expect_equal(
  current = S7::prop(merged_unprefixed, "obs_table")$meta_cell_id,
  target = c(sprintf("left_%02d", 1:6), sprintf("right_%02d", 1:5)),
  info = "mc merge - identifiers stay bare when prefix_ids is FALSE"
)

### dropped columns ------------------------------------------------------------

part_extra <- build_mc(4L, genes_a, seed = 5L)
S7::prop(part_extra$object, "obs_table")[, extra_col := "x"]

expect_warning(
  current = merge_meta_cells(
    inputs = list(a = part_a$object, e = part_extra$object),
    .verbose = FALSE
  ),
  info = "mc merge - warns about observation columns not shared by all inputs"
)

merged_extra <- suppressWarnings(merge_meta_cells(
  inputs = list(a = part_a$object, e = part_extra$object),
  .verbose = FALSE
))

expect_false(
  current = "extra_col" %in% names(S7::prop(merged_extra, "obs_table")),
  info = "mc merge - the non-shared column is dropped"
)

## provenance ------------------------------------------------------------------

expect_true(
  current = S7::prop(merged_int, "is_merged"),
  info = "mc merge - is_merged is set"
)

expect_equal(
  current = S7::prop(merged_int, "other_data")$sources,
  target = c("a", "b"),
  info = "mc merge - the source identifiers land in other_data"
)

assignment_int <- S7::prop(merged_int, "original_assignment")

expect_null(
  current = assignment_int$assignments,
  info = "mc merge - per-cell assignments are dropped, there is no shared space"
)

expect_null(
  current = assignment_int$cells_to_keep,
  info = "mc merge - cells_to_keep is not carried over"
)

expect_equal(
  current = names(assignment_int$per_source),
  target = c("a", "b"),
  info = "mc merge - per_source is keyed by source identifier"
)

expect_equal(
  current = sort(names(assignment_int$per_source$a)),
  target = c("n_cells", "n_metacells", "n_unassigned"),
  info = "mc merge - per_source keeps the scalars only"
)

expect_equal(
  current = assignment_int$n_metacells,
  target = 11L,
  info = "mc merge - n_metacells is the merged count"
)

### shared parent --------------------------------------------------------------

# both hand-built inputs report 100 source cells, so they are treated as sharing
# one obs space and the count is not summed
expect_equal(
  current = assignment_int$n_cells,
  target = 100L,
  info = "mc merge - sources off one parent report that parent's cell count"
)

merged_distinct <- suppressWarnings(merge_meta_cells(
  inputs = list(
    a = part_a$object,
    d = build_mc(4L, genes_a, seed = 6L, n_cells = 120L)$object
  ),
  .verbose = FALSE
))

expect_equal(
  current = S7::prop(merged_distinct, "original_assignment")$n_cells,
  target = 220L,
  info = "mc merge - sources off different parents sum their cell counts"
)

### method -------------------------------------------------------------------

part_other_method <- build_mc(3L, genes_a, seed = 7L)
S7::prop(part_other_method$object, "meta_cell_method") <- "supercells"

expect_error(
  current = merge_meta_cells(
    inputs = list(a = part_a$object, o = part_other_method$object),
    .verbose = FALSE
  ),
  info = "mc merge - inputs from different meta cell methods error"
)

expect_error(
  current = merge_meta_cells(
    inputs = list(part_a$object, part_a$object),
    source_ids = c("same", "same"),
    .verbose = FALSE
  ),
  info = "mc merge - non-unique source identifiers error"
)

## real objects off one parent -------------------------------------------------

fixture <- sc_test_fixture()

sc_object <- sc_test_prepped(
  sc_test_object(test_temp_dir, fixture),
  fixture
)

groups <- c("cell_type_1", "cell_type_2")

mc_inputs <- purrr::map(groups, \(g) {
  subset_obj <- sc_test_prepped(
    SingleCellsSubset(
      sc_object = sc_object,
      grouping_column = "cell_grp",
      group = g
    ),
    fixture
  )
  mc <- generate_bt_meta_cells_sc(
    subset_obj,
    sc_meta_cell_params = params_sc_bt_metacells(target_no_metacells = 25L),
    .verbose = FALSE
  )
  # cached artefacts on the inputs, so the merge can be seen to drop them
  sc_test_prepped(mc, fixture, k = 5L)
})
names(mc_inputs) <- groups

mc_merged <- suppressWarnings(merge_meta_cells(mc_inputs, .verbose = FALSE))

expect_equal(
  current = S7::prop(mc_merged, "dims")[1],
  target = sum(purrr::map_int(mc_inputs, \(x) S7::prop(x, "dims")[1])),
  info = "mc merge - the merged object holds every input meta cell"
)

expect_true(
  current = all(purrr::map_lgl(
    mc_inputs,
    \(x) !is.null(get_pca_factors(x))
  )),
  info = "mc merge - the inputs carry cached artefacts before merging"
)

expect_null(
  current = suppressWarnings(get_pca_factors(mc_merged)),
  info = "mc merge - per-source caches are dropped from the merged object"
)

### refusals -------------------------------------------------------------------

expect_warning(
  current = calc_manifold_metrics(mc_merged),
  info = "mc merge - manifold metrics refuse to run on a merged object"
)

expect_warning(
  current = calc_diffusion_coordinates(
    mc_merged,
    knn_data = get_knn_obj(mc_inputs[[1]]),
    .verbose = FALSE
  ),
  info = "mc merge - diffusion coordinates refuse to run on a merged object"
)

## downstream methods ----------------------------------------------------------

# the point of the merge is that these all survive the trip through Rust

mc_merged <- find_hvg_sc(
  mc_merged,
  hvg_no = 20L,
  .verbose = FALSE
)

expect_true(
  current = length(get_hvg(mc_merged)) == 20L,
  info = "mc merge - HVGs can be identified on the merged object"
)

mc_merged <- calculate_pca_sc(
  mc_merged,
  no_pcs = 5L,
  .verbose = FALSE
)

expect_equal(
  current = nrow(get_pca_factors(mc_merged)),
  target = S7::prop(mc_merged, "dims")[1],
  info = "mc merge - PCA returns one row per merged meta cell"
)

mc_merged <- find_neighbours_sc(
  mc_merged,
  neighbours_params = params_sc_neighbours(knn = list(k = 10L)),
  .verbose = FALSE
)

expect_true(
  current = !is.null(get_knn_obj(mc_merged)),
  info = "mc merge - a kNN graph can be built on the merged object"
)

merged_genes <- S7::prop(mc_merged, "var_table")$gene_id

aucell_res <- aucell_sc(
  mc_merged,
  gs_list = list(
    set_a = merged_genes[1:10],
    set_b = merged_genes[11:20]
  ),
  .verbose = FALSE
)

expect_equal(
  current = dim(aucell_res),
  target = c(S7::prop(mc_merged, "dims")[1], 2L),
  info = "mc merge - AUCell scores one row per merged meta cell"
)

expect_equal(
  current = rownames(aucell_res),
  target = S7::prop(mc_merged, "obs_table")$meta_cell_id,
  info = "mc merge - AUCell rows are labelled with the merged identifiers"
)

nmf_res <- nmf_sc(
  mc_merged,
  k = 3L,
  .verbose = FALSE
)

expect_true(
  current = !is.null(nmf_res),
  info = "mc merge - NMF runs on the merged object"
)

scenic_res <- suppressWarnings(scenic_grn_sc(
  mc_merged,
  tf_ids = merged_genes[1:5],
  .verbose = FALSE
))

expect_true(
  current = !is.null(scenic_res),
  info = "mc merge - SCENIC runs on the merged object"
)

## pipeline entry point --------------------------------------------------------

group_dir <- sc_test_dir("mc_merge_group")

obs_grouped <- data.table::copy(fixture$obs)
obs_grouped[, patient := rep(c("p1", "p2"), length.out = .N)]

sc_grouped <- sc_test_object(group_dir, fixture, obs = obs_grouped)

mc_per_group <- suppressWarnings(meta_cells_per_group(
  sc_grouped,
  group_col = "patient",
  method = "bootstrapped",
  mc_params = list(
    sc_meta_cell_params = params_sc_bt_metacells(target_no_metacells = 25L),
    .verbose = FALSE
  ),
  pipeline = sc_pipeline() %>>%
    step_hvg_sc(hvg_no = fixture$hvg_to_keep, .verbose = FALSE) %>>%
    step_pca_sc(no_pcs = fixture$no_pcs, .verbose = FALSE) %>>%
    step_neighbours_sc(.verbose = FALSE),
  .verbose = FALSE
))

expect_true(
  current = S7::prop(mc_per_group, "is_merged"),
  info = "mc merge - meta_cells_per_group returns a merged object"
)

expect_equal(
  current = sort(unique(S7::prop(mc_per_group, "obs_table")$source_id)),
  target = c("p1", "p2"),
  info = "mc merge - meta_cells_per_group tags each group as a source"
)

sc_test_cleanup(test_temp_dir, group_dir)
