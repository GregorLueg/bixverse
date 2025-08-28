# cor module tests -------------------------------------------------------------

library(magrittr)

## tom -------------------------------------------------------------------------

### data -----------------------------------------------------------------------

cor_data <- c(0.7, 0.3, 0.2, -0.1, 0.25, -0.5)
cor_mat <- rs_upper_triangle_to_dense(cor_data, 1L, 4)

tom_signed_v1 <- c(
  0.306382979,
  0.05,
  0.081818182,
  -0.005357143,
  0.162962963,
  -0.19375
)

tom_signed_v2 <- c(
  0.3536364,
  0.1113636,
  0.1058140,
  -0.02875,
  0.1681818,
  -0.2427083
)

tom_unsigned_v1 <- c(
  0.3319149,
  0.1807692,
  0.1909091,
  0.1553571,
  0.1629630,
  0.2437500
)

tom_unsigned_v2 <- c(
  0.3645455,
  0.1886364,
  0.1755814,
  0.1337500,
  0.1681818,
  0.2677083
)

### tests ----------------------------------------------------------------------

tom_v1_signed_mat <- calculate_tom(
  cor_mat = cor_mat,
  signed = TRUE,
  version = "v1"
)

tom_v2_signed_mat <- calculate_tom(
  cor_mat = cor_mat,
  signed = TRUE,
  version = "v2"
)

tom_v1_unsigned_mat <- calculate_tom(
  cor_mat = cor_mat,
  signed = FALSE,
  version = "v1"
)

tom_v2_unsigned_mat <- calculate_tom(
  cor_mat = cor_mat,
  signed = FALSE,
  version = "v2"
)

expect_equal(
  current = rs_dense_to_upper_triangle(tom_v1_signed_mat, 1),
  target = tom_signed_v1,
  info = paste(
    "TOM calculation version 1 - signed"
  ),
  tolerance = 10e-6
)

expect_equal(
  current = rs_dense_to_upper_triangle(tom_v2_signed_mat, 1),
  target = tom_signed_v2,
  info = paste(
    "TOM calculation version 2 - signed"
  ),
  tolerance = 10e-6
)

expect_equal(
  current = rs_dense_to_upper_triangle(tom_v1_unsigned_mat, 1),
  target = tom_unsigned_v1,
  info = paste(
    "TOM calculation version 1 - unsigned"
  ),
  tolerance = 10e-6
)

expect_equal(
  current = rs_dense_to_upper_triangle(tom_v2_unsigned_mat, 1),
  target = tom_unsigned_v2,
  info = paste(
    "TOM calculation version 2 - unsigned"
  ),
  tolerance = 10e-6
)

## coremo ----------------------------------------------------------------------

### synthetic data -------------------------------------------------------------

test_data <- rs_generate_bulk_rnaseq(
  num_samples = 100L,
  num_genes = 1000L,
  seed = 123L,
  add_modules = TRUE,
  module_sizes = c(100L, 100L, 100L)
)

norm_counts <- edgeR::cpm(test_data$counts, log = TRUE)

rownames(norm_counts) <- sprintf("gene_%i", 1:1000)
colnames(norm_counts) <- sprintf("sample_%i", 1:100)

data <- t(norm_counts)
data_bad <- data
rownames(data_bad) <- NULL
colnames(data_bad) <- NULL

meta_data <- data.table::data.table(
  sample_id = rownames(data)
)

meta_data_bad <- data.table::data.table(
  sample_id = sprintf("sample_%i", 101:200)
)

### expected data --------------------------------------------------------------

expected_coremo_modules <- qs2::qs_read("./test_data/coremo_modules.qs")

## tests -----------------------------------------------------------------------

### object generation ----------------------------------------------------------

expect_error(
  current = bulk_coexp(raw_data = data_bad, meta_data = meta_data),
  info = "bulk_coexp - bad matrix"
)

expect_error(
  current = bulk_coexp(raw_data = data, meta_data = meta_data_bad),
  info = "bulk_coexp - bad metadata"
)

cor_test <- bulk_coexp(raw_data = data, meta_data = meta_data)

expect_true(
  S7::S7_inherits(cor_test, bulk_coexp),
  info = "bulk_coexp - correct main class"
)
expect_true(
  S7::S7_inherits(cor_test, bixverse_base_class),
  info = "bulk_coexp - correct base class"
)

### hvg ------------------------------------------------------------------------

expect_warning(
  current = plot_hvgs(cor_test),
  info = "bulk_coexp - warning without hvg genes"
)

cor_test <- preprocess_bulk_coexp(
  cor_test,
  mad_threshold = 0.4,
  .verbose = FALSE
)

expect_true(
  current = sum(cor_test@processed_data$feature_meta$hvg) == 737,
  info = "bulk_coexp - correct number of HVG genes with mad thresholding"
)

hvg_plot <- plot_hvgs(cor_test)

expect_true(
  "ggplot" %in% class(hvg_plot),
  info = "bulk_coexp - plotting of hvg working"
)

cor_test <- preprocess_bulk_coexp(cor_test, hvg = 0.5, .verbose = FALSE)

expect_true(
  current = sum(cor_test@processed_data$feature_meta$hvg) == 500,
  info = "bulk_coexp - correct number of HVG genes with proportion"
)

### correlations and epsilon search --------------------------------------------

cor_test <- cor_module_processing(
  cor_test,
  cor_method = "spearman",
  .verbose = FALSE
)

cor_matrix <- S7::prop(cor_test, "processed_data")[[
  "correlation_res"
]]$get_sym_matrix(
  .verbose = FALSE
)

expect_equal(
  current = dim(cor_matrix),
  target = c(500, 500),
  info = paste(
    "CoReMo - cor matrix dimensions"
  )
)

expect_true(
  is.null(cor_test@outputs$epsilon_data),
  info = "CoReMo - no epsilon data found"
)

expect_warning(
  current = plot_epsilon_res(cor_test),
  info = "CoReMo - warning without plot epsilon"
)

cor_test <- cor_module_check_epsilon(
  cor_test,
  rbf_func = "gaussian",
  .verbose = FALSE
)

expect_true(
  !is.null(cor_test@outputs$epsilon_data),
  info = "CoReMo - epsilon data generated"
)

epsilon_plot <- plot_epsilon_res(cor_test)

expect_true(
  "ggplot" %in% class(epsilon_plot),
  info = "CoReMo - plotting of epsilon results working"
)

chosen_epsilon <- cor_test@outputs$epsilon_data[
  r2_vals == max(r2_vals),
  epsilon
]

expect_equal(
  current = chosen_epsilon,
  target = 8,
  info = "CoReMo - correct epsilon"
)

### clustering and stability ---------------------------------------------------

expect_warning(
  current = plot_optimal_cuts(cor_test),
  info = "CoReMo - warning without plot optimal cuts"
)

cor_test <- cor_module_coremo_clustering(cor_test, .verbose = FALSE)

expect_true(
  !is.null(cor_test@outputs$optimal_cuts),
  info = "CoReMo - optimal cuts data generated"
)

plot_cuts <- plot_optimal_cuts(cor_test)

expect_true(
  "ggplot" %in% class(plot_cuts),
  info = "CoReMo - plotting of optimal cuts working"
)

cor_test <- cor_module_coremo_stability(cor_test, .verbose = FALSE)

final_module_data <- cor_test@outputs$final_modules

expect_equal(
  current = final_module_data,
  target = expected_coremo_modules,
  info = "CoReMo - final module data"
)

## graph-based clustering ------------------------------------------------------

### expected data --------------------------------------------------------------

expected_cor_graph_res <- qs2::qs_read(
  "./test_data/cor_graph_final_res.qs"
)

### regenerate the data --------------------------------------------------------

cor_test <- bulk_coexp(raw_data = data, meta_data = meta_data) %>%
  preprocess_bulk_coexp(., hvg = 0.5, .verbose = FALSE) %>%
  cor_module_processing(., cor_method = "spearman", .verbose = FALSE)

### epsilon results ------------------------------------------------------------

cor_test <- cor_module_check_epsilon(
  cor_test,
  rbf_func = "bump",
  .verbose = FALSE
)

epsilon_results <- get_epsilon_res(cor_test)

expect_true(
  current = checkmate::test_data_table(epsilon_results),
  info = "graph-based clustering - epsilon data getter"
)

chosen_pump_epsilon <- epsilon_results[r2_vals == max(r2_vals), epsilon]

expect_equal(
  current = chosen_pump_epsilon,
  target = 10,
  info = "graph-based clustering - correct epsilon"
)

### resolution iterations ------------------------------------------------------

expect_warning(
  current = get_resolution_res(cor_test),
  info = paste(
    "Testing the warning of missing resolution results."
  )
)

cor_test <- cor_module_graph_check_res(
  object = cor_test,
  .verbose = FALSE
)

resolution_results <- get_resolution_res(cor_test)

expect_true(
  current = checkmate::test_data_table(resolution_results),
  info = "graph-based clustering - resolution data getter"
)

### finalisation of the clusters -----------------------------------------------

cor_test <- cor_module_graph_final_modules(
  object = cor_test,
  .verbose = FALSE
)

final_cor_graph_res <- get_results(cor_test)

expect_equivalent(
  current = final_cor_graph_res,
  target = expected_cor_graph_res,
  info = paste(
    "Testing expected final results from the graph-based cor results"
  )
)
