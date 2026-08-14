# sc cache provenance state ----------------------------------------------------

source("helper_sc.R", local = TRUE)

set.seed(123L)

test_temp_dir <- sc_test_dir("cache_state")

## fixture ---------------------------------------------------------------------

single_cell_test_data <- sc_test_fixture()

sc_qc_param <- sc_test_qc_params(single_cell_test_data, target_size = 1000)

sc_object <- sc_test_object(
  test_temp_dir,
  single_cell_test_data,
  sc_qc_param = sc_qc_param
)

sc_object <- find_hvg_sc(sc_object, hvg_no = 30L, .verbose = FALSE)
sc_object <- calculate_pca_sc(sc_object, no_pcs = 10L, .verbose = FALSE)
sc_object <- umap_sc(sc_object, use_knn = FALSE, .verbose = FALSE)
sc_object <- find_neighbours_sc(sc_object, .verbose = FALSE)

## happy path ------------------------------------------------------------------

status <- get_sc_cache_status(sc_object)

expect_true(
  current = data.table::is.data.table(status),
  info = "cache status returns a data.table"
)

expect_true(
  current = all(
    c("modality", "artefact", "name", "stamped", "stale", "reason", "id") %in%
      names(status)
  ),
  info = "cache status has the documented columns"
)

expect_true(
  current = nrow(status) == 4L,
  info = "cache status covers pca, umap, knn and snn"
)

expect_true(
  current = all(status$stamped),
  info = "freshly computed artefacts are all stamped"
)

expect_false(
  current = any(status$stale),
  info = "a clean pipeline leaves nothing stale"
)

expect_true(
  current = !anyDuplicated(status$id) && all(nchar(status$id) == 16L),
  info = "stamp ids are unique and 16 characters"
)

# the row index is computed outside the `[` on purpose: `artefact` and `name`
# are also column names, and data.table would resolve them to the columns
.stamp_row <- function(status, artefact, name = NA_character_) {
  idx <- if (!is.na(name)) {
    which(status$name == name)
  } else {
    which(status$artefact == artefact)
  }
  status[idx, ]
}

pca_id <- .stamp_row(status, "pca")$id
knn_id <- .stamp_row(status, "knn")$id

expect_equal(
  current = .stamp_row(status, "embedding", "umap")$from[[1]],
  target = pca_id,
  info = "the umap records the pca it was built on"
)

expect_equal(
  current = .stamp_row(status, "knn")$from[[1]],
  target = pca_id,
  info = "the knn records the pca it was built on"
)

expect_equal(
  current = .stamp_row(status, "snn")$from[[1]],
  target = knn_id,
  info = "the snn records the knn it was built on"
)

## attribute hygiene -----------------------------------------------------------

expect_null(
  current = attr(get_pca_factors(sc_object), "bixverse_stamp"),
  info = "the pca getter does not leak the stamp to user code"
)

expect_null(
  current = attr(get_embedding(sc_object, "umap"), "bixverse_stamp"),
  info = "the embedding getter does not leak the stamp"
)

expect_null(
  current = attr(get_knn_obj(sc_object), "bixverse_stamp"),
  info = "the knn getter does not leak the stamp"
)

expect_null(
  current = attr(get_snn_graph(sc_object), "bixverse_stamp"),
  info = "the snn getter does not leak the stamp"
)

## the reported foot gun -------------------------------------------------------

all_cells <- get_cell_names(sc_object, filtered = TRUE)
fewer_cells <- all_cells[-(1:5)]

expect_warning(
  current = {
    sc_filtered <- set_cells_to_keep(sc_object, fewer_cells)
  },
  pattern = "stale",
  info = "filtering cells warns about the artefacts it just invalidated"
)

filtered_status <- get_sc_cache_status(sc_filtered)

expect_true(
  current = all(filtered_status$stale),
  info = "every artefact is stale once the cell set moved"
)

expect_warning(
  current = get_pca_factors(sc_filtered),
  pattern = "stale",
  info = "the pca getter warns on a stale artefact"
)

expect_error(
  current = find_neighbours_sc(sc_filtered, .verbose = FALSE),
  pattern = "stale",
  info = "the hard tier errors rather than building on a stale embedding"
)

# re-run the pca only. this is the exact scenario from the bug report: the
# cell set and the pca now agree, but the knn was never regenerated.
sc_repaired <- calculate_pca_sc(sc_filtered, no_pcs = 10L, .verbose = FALSE)

repaired_status <- get_sc_cache_status(sc_repaired)

expect_false(
  current = .stamp_row(repaired_status, "pca")$stale,
  info = "re-running the pca makes it fresh again"
)

expect_true(
  current = .stamp_row(repaired_status, "knn")$stale,
  info = "the knn stays stale after the pca alone was re-run"
)

expect_true(
  current = .stamp_row(repaired_status, "snn")$stale,
  info = "the snn is stale through its upstream knn, not directly"
)

expect_equal(
  current = .stamp_row(repaired_status, "snn")$reason,
  target = "the cell set changed since it was computed",
  info = "the snn fails on its own cell hash before the chain is consulted"
)

expect_error(
  current = run_paga_sc(
    sc_repaired,
    cluster_col = "cell_grp",
    .verbose = FALSE
  ),
  pattern = "stale",
  info = "trajectory refuses to run on a stale knn"
)

sc_repaired <- find_neighbours_sc(sc_repaired, .verbose = FALSE)

final_status <- get_sc_cache_status(sc_repaired)

expect_false(
  current = any(final_status$stale[final_status$artefact != "embedding"]),
  info = "re-running the neighbours clears the chain"
)

expect_equal(
  current = .stamp_row(final_status, "snn")$from[[1]],
  target = .stamp_row(final_status, "knn")$id,
  info = "the snn points at the new knn id, not the old one"
)

## same cells, different pca ---------------------------------------------------

# the case a cell-set hash on its own cannot catch: nothing was filtered, the
# pca was simply re-run with a different number of components, and the kNN
# built on the old one was never regenerated
sc_repca <- calculate_pca_sc(sc_object, no_pcs = 5L, .verbose = FALSE)

repca_status <- get_sc_cache_status(sc_repca)

expect_false(
  current = .stamp_row(repca_status, "pca")$stale,
  info = "the re-run pca is fresh"
)

expect_equal(
  current = .stamp_row(repca_status, "knn")$reason,
  target = "the artefact it was derived from was re-computed or removed",
  info = "the knn is stale because its parent pca minted a new id"
)

expect_equal(
  current = .stamp_row(repca_status, "snn")$reason,
  target = "its upstream `rna:knn` is stale",
  info = "the snn is stale transitively, through an otherwise valid knn"
)

expect_error(
  current = find_clusters_sc(sc_repca),
  pattern = "stale",
  info = "clustering refuses a graph whose provenance chain is broken"
)

## severity option -------------------------------------------------------------

old_opts <- options(bixverse.cache_check = "error")

expect_error(
  current = get_pca_factors(sc_filtered),
  pattern = "stale",
  info = "the error mode promotes the getter warning"
)

options(bixverse.cache_check = "none")

# the knn getter, not the pca one: with the checks off the stale pca still
# blows up on its own rownames assignment, which is exactly what "none" buys
# you. silencing the stamp machinery does not make a mis-shaped matrix safe.
expect_silent(
  current = get_knn_obj(sc_filtered),
  info = "the none mode silences the getter"
)

expect_silent(
  current = assert_sc_state(sc_filtered, artefacts = "pca"),
  info = "the none mode also disables the hard tier"
)

options(old_opts)

## legacy objects --------------------------------------------------------------

# an object built before provenance stamping existed carries no attributes
sc_legacy <- sc_object
legacy_cache <- S7::prop(sc_legacy, "sc_cache")
legacy_cache$pca_factors <- bixverse:::.drop_stamp(legacy_cache$pca_factors)
legacy_cache$knn <- bixverse:::.drop_stamp(legacy_cache$knn)
legacy_cache$snn_graph <- bixverse:::.drop_stamp(legacy_cache$snn_graph)
legacy_cache$other_embeddings$umap <- bixverse:::.drop_stamp(
  legacy_cache$other_embeddings$umap
)
S7::prop(sc_legacy, "sc_cache") <- legacy_cache

legacy_status <- get_sc_cache_status(sc_legacy)

expect_false(
  current = any(legacy_status$stamped),
  info = "stripped artefacts report as unstamped"
)

expect_false(
  current = any(legacy_status$stale),
  info = "unknown provenance is not the same as known bad"
)

expect_silent(
  current = get_pca_factors(sc_legacy),
  info = "legacy artefacts still read without complaint"
)

# and they stay silent even when the cell set genuinely moved, because
# stamp-only means there is nothing to compare against
sc_legacy_filtered <- suppressWarnings(
  set_cells_to_keep(sc_legacy, fewer_cells)
)

expect_silent(
  current = get_knn_obj(sc_legacy_filtered),
  info = "unstamped artefacts pass even after the cell set moved"
)

## removal ---------------------------------------------------------------------

sc_removed <- remove_knn(sc_object)

removed_status <- get_sc_cache_status(sc_removed)

expect_true(
  current = !"knn" %in% removed_status$artefact,
  info = "removing the knn removes its stamp with it"
)

expect_true(
  current = .stamp_row(removed_status, "snn")$stale,
  info = "the snn goes stale once the knn it derives from is gone"
)

## round trip ------------------------------------------------------------------

save_sc_exp_to_disk(sc_object, type = "rds")

sc_reloaded <- load_existing(SingleCells(dir_data = test_temp_dir))

expect_equal(
  current = get_sc_cache_status(sc_reloaded)$id,
  target = status$id,
  info = "stamp ids survive a save and load round trip"
)

## print -----------------------------------------------------------------------

expect_stdout(
  current = print(sc_object),
  pattern = "Stale artefacts: none",
  info = "print reports a clean object as clean"
)

expect_stdout(
  current = print(sc_filtered),
  pattern = "Stale artefacts: rna:pca",
  info = "print names the stale artefacts"
)

## reset -----------------------------------------------------------------------

# tinytest runs non-interactively, which is the branch that matters: there is
# nobody to confirm to, so a destructive reset must refuse rather than guess
if (!interactive()) {
  expect_error(
    current = reset_cells_to_keep(sc_filtered),
    pattern = "force",
    info = "a non-interactive reset without force errors instead of prompting"
  )
}

sc_reset <- reset_cells_to_keep(sc_filtered, force = TRUE)

expect_equal(
  current = length(get_cells_to_keep(sc_reset)),
  target = as.integer(get_sc_rust_ptr(sc_reset)$get_shape()[1]),
  info = "reset restores every cell in the binary"
)

expect_equal(
  current = length(get_sc_duckdb(sc_reset)$get_cells_to_keep()),
  target = as.integer(get_sc_rust_ptr(sc_reset)$get_shape()[1]),
  info = "the duckdb agrees with the map after a reset"
)

expect_true(
  current = nrow(get_sc_cache_status(sc_reset)) == 0L,
  info = "reset leaves an empty cache"
)

expect_warning(
  current = get_hvg(sc_reset),
  info = "reset clears the hvg selection, which is cell set dependent"
)

# a reset object must be genuinely usable, not half wired
sc_reset <- find_hvg_sc(sc_reset, hvg_no = 30L, .verbose = FALSE)
sc_reset <- calculate_pca_sc(sc_reset, no_pcs = 10L, .verbose = FALSE)
sc_reset <- find_neighbours_sc(sc_reset, .verbose = FALSE)

expect_false(
  current = any(get_sc_cache_status(sc_reset)$stale),
  info = "the full pipeline runs clean on a reset object"
)

# nothing cached and nothing to lose: no prompt, no error, even without force
sc_empty <- reset_cells_to_keep(
  reset_cells_to_keep(sc_reset, force = TRUE)
)

expect_true(
  current = nrow(get_sc_cache_status(sc_empty)) == 0L,
  info = "resetting an already pristine object is a silent no-op"
)

## subset and meta cells -------------------------------------------------------

sc_subset <- SingleCellsSubset(
  sc_object,
  "cell_grp",
  unlist(sc_object[["cell_grp"]])[1]
)

expect_error(
  current = reset_cells_to_keep(sc_subset, force = TRUE),
  pattern = "subset",
  info = "a subset refuses to be reset out of being a subset"
)

expect_true(
  current = nrow(get_sc_cache_status(sc_subset)) == 0L,
  info = "a fresh subset starts with an empty cache"
)

# clean up ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir)
