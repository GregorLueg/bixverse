# sc bonsai --------------------------------------------------------------------

source("helper_sc.R", local = TRUE)

test_temp_dir <- sc_test_dir("bonsai")

## fixtures --------------------------------------------------------------------

fixture <- sc_test_fixture()
sc_object <- sc_test_object(test_temp_dir, fixture)

cell_names <- get_cell_names(sc_object, filtered = TRUE)
all_genes <- get_gene_names(sc_object)

tree <- bonsai_sc(sc_object, .verbose = FALSE)

# tests ------------------------------------------------------------------------

## tree structure --------------------------------------------------------------

expect_true(
  inherits(tree, "BonsaiTree"),
  info = "bonsai: returns a BonsaiTree"
)

expect_equal(
  tree$nodes[(is_leaf)]$cell_id,
  cell_names,
  info = "bonsai: one leaf per cell, in the object's cell order"
)

expect_equal(
  sum(is.na(tree$nodes$parent)),
  1L,
  info = "bonsai: exactly one root"
)

expect_true(
  all(tree$nodes[!is.na(parent)]$parent > tree$nodes[!is.na(parent)]$node),
  info = "bonsai: parents come after their children"
)

expect_true(
  all(is.finite(c(tree$nodes$x, tree$nodes$y))),
  info = "bonsai: finite coordinates"
)

expect_true(
  length(tree$genes_used) > 0,
  info = "bonsai: all genes in, some pass the signal-to-noise filter"
)

expect_error(
  bonsai_sc(
    sc_object,
    bonsai_params = params_sc_bonsai(min_signal_to_noise = 1e12),
    .verbose = FALSE
  ),
  pattern = "retained|feature",
  info = "bonsai: a threshold no gene passes errors cleanly"
)

expect_equal(
  sort(c(tree$genes_used, tree$genes_dropped)),
  sort(all_genes),
  info = "bonsai: used and dropped genes partition all genes"
)

expect_equal(
  tree$timings$stage,
  c("sanity", "ingest", "bonsai", "layout", "total"),
  info = "bonsai: every stage timed, plus the total"
)

expect_true(
  all(tree$timings$seconds >= 0) &&
    tree$timings[stage == "total", seconds] >=
      sum(tree$timings[stage != "total", seconds]),
  info = "bonsai: the total covers the stages"
)

expect_stdout(
  print(tree),
  "BonsaiTree:",
  info = "bonsai: print method registered"
)

## relayout --------------------------------------------------------------------

tree_hyp <- relayout_bonsai(tree, layout = "dendrogram", hyperbolic = TRUE)

expect_equal(
  tree_hyp$nodes[(is_leaf)]$cell_id,
  cell_names,
  info = "bonsai relayout: leaves keep their cells"
)

expect_equal(
  nrow(tree_hyp$nodes),
  nrow(tree$nodes),
  info = "bonsai relayout: same number of nodes"
)

expect_true(
  all(tree_hyp$nodes$x^2 + tree_hyp$nodes$y^2 <= 1),
  info = "bonsai relayout: hyperbolic layout stays in the unit disk"
)

expect_true(
  tree_hyp$layout == "dendrogram" && tree_hyp$hyperbolic,
  info = "bonsai relayout: layout recorded"
)

## plotting --------------------------------------------------------------------

expect_true(
  inherits(
    plot(
      tree,
      colour_by = get_sc_obs(sc_object, filtered = TRUE)$cell_grp
    ),
    "ggplot"
  ),
  info = "bonsai plot: returns a ggplot"
)

expect_true(
  inherits(plot(tree, layout = "dendrogram"), "ggplot"),
  info = "bonsai plot: relayouts on the fly"
)

expect_error(
  plot(tree, colour_by = 1:3),
  info = "bonsai plot: colour_by has to cover every cell"
)

## embedding -------------------------------------------------------------------

sc_object <- set_bonsai_embedding(sc_object, tree)
embd <- get_embedding(sc_object, "bonsai")

expect_equal(
  dim(embd),
  c(length(cell_names), 2L),
  info = "bonsai embedding: cells x 2"
)

expect_equal(
  unname(embd[, 1]),
  tree$nodes[(is_leaf)]$x,
  info = "bonsai embedding: leaf x coordinates"
)

## parameters and edge cases ---------------------------------------------------

expect_error(
  params_sc_bonsai(layout = "circle"),
  info = "bonsai params: unknown layout"
)

expect_error(
  bonsai_sc(
    sc_object,
    bonsai_params = params_sc_bonsai(variance_rule = "fixed"),
    .verbose = FALSE
  ),
  info = "bonsai: fixed variance rule without a variance errors in Rust"
)

tree_fixed <- bonsai_sc(
  sc_object,
  bonsai_params = params_sc_bonsai(
    variance_rule = "fixed",
    fixed_variance = 1.0,
    layout = "dendrogram"
  ),
  .verbose = FALSE
)

expect_equal(
  tree_fixed$layout,
  "dendrogram",
  info = "bonsai: layout from the params"
)

expect_equal(
  bonsai_sc(sc_object, .verbose = FALSE)$loglik,
  tree$loglik,
  info = "bonsai: same seed, same tree"
)

sc_object <- find_hvg_sc(
  sc_object,
  hvg_no = fixture$hvg_to_keep,
  .verbose = FALSE
)
hvg_names <- unname(get_gene_names_from_idx(sc_object, get_hvg(sc_object)))
tree_hvg <- bonsai_sc(
  sc_object,
  hvg = get_hvg(sc_object) + 1L,
  .verbose = FALSE
)

expect_true(
  length(tree_hvg$genes_used) > 0 && all(tree_hvg$genes_used %in% hvg_names),
  info = "bonsai hvg: genes used come from the given genes"
)

expect_equal(
  sort(c(tree_hvg$genes_used, tree_hvg$genes_dropped)),
  sort(hvg_names),
  info = "bonsai hvg: used and dropped genes partition the given genes"
)

## metacells -------------------------------------------------------------------

mc_dir <- sc_test_dir("bonsai_mc")
mc_source <- sc_test_prepped(sc_test_object(mc_dir, fixture), fixture)
mc_object <- generate_bt_meta_cells_sc(
  mc_source,
  sc_meta_cell_params = params_sc_bt_metacells(target_no_metacells = 200L),
  .verbose = FALSE
)
mc_obs <- mc_object[[]]

tree_mc <- bonsai_sc(mc_object, .verbose = FALSE)

expect_equal(
  tree_mc$nodes[(is_leaf)]$cell_id,
  mc_obs$meta_cell_id,
  info = "bonsai mc: one leaf per metacell, in the object's order"
)

expect_equal(
  tree_mc$nodes[(is_leaf)]$n_cells,
  as.integer(mc_obs$no_originating_cells),
  info = "bonsai mc: leaves carry their metacell sizes"
)

expect_true(
  all(is.na(tree_mc$nodes[!(is_leaf)]$n_cells)) &&
    all(tree$nodes[(is_leaf)]$n_cells == 1L),
  info = "bonsai: ancestors have no size, single cells have one"
)

expect_equal(
  sum(is.na(tree_mc$nodes$parent)),
  1L,
  info = "bonsai mc: exactly one root"
)

expect_stdout(
  print(tree_mc),
  "cells\\)",
  info = "bonsai mc: print reports the cells behind the leaves"
)

expect_true(
  inherits(plot(tree_mc, size_by_cells = TRUE), "ggplot"),
  info = "bonsai mc: plot sized by metacell size"
)

mc_object <- set_bonsai_embedding(mc_object, tree_mc)

expect_equal(
  dim(get_embedding(mc_object, "bonsai")),
  c(nrow(mc_obs), 2L),
  info = "bonsai mc: embedding is metacells x 2"
)

# clean up ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir, mc_dir)
