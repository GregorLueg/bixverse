# sc processing ----------------------------------------------------------------

source("helper_sc.R", local = TRUE)

set.seed(123L)

test_temp_dir <- sc_test_dir("plot_helpers")

## testing parameters ----------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
# hvg
hvg_to_keep <- 30L
# pca
no_pcs <- 15L

## synthetic test data ---------------------------------------------------------

single_cell_test_data <- sc_test_fixture(
  min_lib_size = min_lib_size,
  min_genes_exp = min_genes_exp,
  min_cells_exp = min_cells_exp,
  hvg_to_keep = hvg_to_keep,
  no_pcs = no_pcs
)

genes_pass <- single_cell_test_data$genes_pass
cells_pass <- single_cell_test_data$cells_pass

counts_filtered <- single_cell_test_data$counts[cells_pass, genes_pass]

sc_qc_param <- sc_test_qc_params(single_cell_test_data, target_size = 1000)

## underlying class ------------------------------------------------------------

sc_object <- sc_test_object(test_temp_dir, single_cell_test_data)

# tests ------------------------------------------------------------------------

## extract_dot_plot_data -------------------------------------------------------

test_features <- sprintf("gene_%03d", 1:5)

dot_dt <- extract_dot_plot_data(
  object = sc_object,
  features = test_features,
  grouping_variable = "cell_grp",
  scale_exp = TRUE
)

n_clusters <- length(unique(
  unlist(sc_object[["cell_grp"]], use.names = FALSE)
))

expect_true(
  current = data.table::is.data.table(dot_dt),
  info = "extract_dot_plot_data returns a data.table"
)

expect_equal(
  current = nrow(dot_dt),
  target = length(test_features) * n_clusters,
  info = "extract_dot_plot_data returns correct number of rows"
)

expect_true(
  current = all(
    c("gene", "group", "mean_exp", "scaled_exp", "pct_exp") %in%
      names(dot_dt)
  ),
  info = "extract_dot_plot_data has expected columns"
)

expect_true(
  current = is.factor(dot_dt$gene) && is.factor(dot_dt$group),
  info = "extract_dot_plot_data gene and group are factors"
)

expect_true(
  current = all(dot_dt$pct_exp >= 0 & dot_dt$pct_exp <= 100),
  info = "extract_dot_plot_data pct_exp in [0, 100]"
)

expect_true(
  current = all(dot_dt$scaled_exp >= 0 & dot_dt$scaled_exp <= 1),
  info = "extract_dot_plot_data scaled_exp in [0, 1] after min-max scaling"
)

# check that per gene at least one value is 1 (the max)
expect_true(
  current = all(
    dot_dt[, max(scaled_exp), by = gene][V1 > 0, V1] == 1
  ),
  info = "extract_dot_plot_data scaled_exp max is 1 per gene (where non-zero)"
)

# unscaled version
dot_dt_unscaled <- extract_dot_plot_data(
  object = sc_object,
  features = test_features,
  grouping_variable = "cell_grp",
  scale_exp = FALSE
)

expect_equal(
  current = dot_dt_unscaled$scaled_exp,
  target = dot_dt_unscaled$mean_exp,
  info = "extract_dot_plot_data without scaling scaled_exp equals mean_exp"
)

## extract_gene_expression -----------------------------------------------------

expr_dt <- extract_gene_expression(
  object = sc_object,
  features = test_features
)

expect_true(
  current = data.table::is.data.table(expr_dt),
  info = "extract_gene_expression returns a data.table"
)

expect_equal(
  current = nrow(expr_dt),
  target = length(get_cell_names(sc_object)),
  info = "extract_gene_expression has one row per cell"
)

expect_true(
  current = all(c("cell_id", test_features) %in% names(expr_dt)),
  info = "extract_gene_expression has cell_id and gene columns"
)

expect_true(
  current = all(sapply(test_features, function(g) is.numeric(expr_dt[[g]]))),
  info = "extract_gene_expression gene columns are numeric"
)

# with obs columns
expr_dt_obs <- extract_gene_expression(
  object = sc_object,
  features = test_features,
  obs_cols = c("batch_index", "cell_grp")
)

expect_true(
  current = all(c("batch_index", "cell_grp") %in% names(expr_dt_obs)),
  info = "extract_gene_expression obs_cols are appended"
)

expect_equal(
  current = expr_dt_obs$cell_grp,
  target = unlist(sc_object[["cell_grp"]], use.names = FALSE),
  info = "extract_gene_expression obs_cols values match object"
)

# with scaling
expr_dt_scaled <- extract_gene_expression(
  object = sc_object,
  features = test_features,
  scale = TRUE
)

expect_true(
  current = all(sapply(test_features, function(g) {
    abs(mean(expr_dt_scaled[[g]])) < 0.1
  })),
  info = "extract_gene_expression scaled values have near-zero mean"
)

# with scaling and clipping
expr_dt_clipped <- extract_gene_expression(
  object = sc_object,
  features = test_features,
  scale = TRUE,
  clip = 2.0
)

expect_true(
  current = all(sapply(test_features, function(g) {
    all(expr_dt_clipped[[g]] >= -2.0 & expr_dt_clipped[[g]] <= 2.0)
  })),
  info = "extract_gene_expression clipped values within [-clip, clip]"
)

## unmatched features ----------------------------------------------------------

# `get_gene_indices()` drops what it cannot match. Anything labelling a Rust
# result by position has to follow that, or the labels slide off the values.

mixed_features <- c(test_features[1:2], "not_a_gene", test_features[3])

expect_warning(
  current = extract_gene_expression(
    object = sc_object,
    features = mixed_features
  ),
  info = "extract_gene_expression warns about features it cannot match"
)

expr_dt_mixed <- suppressWarnings(extract_gene_expression(
  object = sc_object,
  features = mixed_features
))

expect_equal(
  current = setdiff(names(expr_dt_mixed), "cell_id"),
  target = test_features[c(1, 2, 3)],
  info = "extract_gene_expression only names the columns it actually got"
)

expr_dt_clean <- extract_gene_expression(
  object = sc_object,
  features = test_features[c(1, 2, 3)]
)

expect_equal(
  current = expr_dt_mixed,
  target = expr_dt_clean,
  info = "extract_gene_expression drops the unmatched gene without shifting"
)

dot_dt_mixed <- suppressWarnings(extract_dot_plot_data(
  object = sc_object,
  features = mixed_features,
  grouping_variable = "cell_grp",
  scale_exp = TRUE
))

expect_equal(
  current = levels(dot_dt_mixed$gene),
  target = test_features[c(1, 2, 3)],
  info = "extract_dot_plot_data only labels the genes it actually got"
)

dot_dt_clean <- extract_dot_plot_data(
  object = sc_object,
  features = test_features[c(1, 2, 3)],
  grouping_variable = "cell_grp",
  scale_exp = TRUE
)

expect_equal(
  current = dot_dt_mixed,
  target = dot_dt_clean,
  info = "extract_dot_plot_data drops the unmatched gene without shifting"
)

# two dropped out of four used to recycle cleanly and mislabel in silence
two_missing <- c(test_features[1], "nope_1", test_features[2], "nope_2")

dot_dt_two <- suppressWarnings(extract_dot_plot_data(
  object = sc_object,
  features = two_missing,
  grouping_variable = "cell_grp",
  scale_exp = TRUE
))

expect_equal(
  current = levels(dot_dt_two$gene),
  target = test_features[1:2],
  info = "extract_dot_plot_data survives an evenly recycling drop"
)

expect_equal(
  current = nrow(dot_dt_two),
  target = 2L * n_clusters,
  info = "extract_dot_plot_data returns one row per surviving gene and group"
)

# clean up ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir)
