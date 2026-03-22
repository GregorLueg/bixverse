# sc processing ----------------------------------------------------------------

set.seed(123L)

test_temp_dir <- file.path(
  tempdir(),
  paste0("test_", format(Sys.time(), "%Y%m%d_%H%M%S_"), sample(1000:9999, 1))
)
dir.create(test_temp_dir, recursive = TRUE)

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

single_cell_test_data <- generate_single_cell_test_data()

genes_pass <- which(
  Matrix::colSums(single_cell_test_data$counts != 0) >= min_cells_exp
)

cells_pass <- which(
  (Matrix::rowSums(single_cell_test_data$counts[, genes_pass]) >=
    min_lib_size) &
    (Matrix::rowSums(single_cell_test_data$counts[, genes_pass] != 0) >=
      min_genes_exp)
)

expect_true(
  current = length(genes_pass) > 80 & length(genes_pass) != 100,
  info = "sc processing - sensible amount of genes pass"
)

expect_true(
  current = length(cells_pass) > 800 & length(cells_pass) != 1000,
  info = "sc processing - sensible amount of cells pass"
)

counts_filtered <- single_cell_test_data$counts[cells_pass, genes_pass]

sc_qc_param <- params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp,
  target_size = 1000
)

## underlying class ------------------------------------------------------------

sc_object <- SingleCells(dir_data = test_temp_dir)

sc_object <- load_r_data(
  object = sc_object,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  var = single_cell_test_data$var,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = min_genes_exp,
    min_lib_size = min_lib_size,
    min_cells = min_cells_exp
  ),
  streaming = FALSE,
  .verbose = FALSE
)

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

expect_equal(
  current = levels(dot_dt$gene),
  target = rev(test_features),
  info = "extract_dot_plot_data gene factor levels in reverse input order"
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
