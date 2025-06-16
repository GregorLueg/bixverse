# cor module tests -------------------------------------------------------------

library(magrittr)

## data ------------------------------------------------------------------------

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

### tom ------------------------------------------------------------------------

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

### synthetic data -------------------------------------------------------------

syn_data <- synthetic_signal_matrix()

data <- t(syn_data$mat)

meta_data <- data.table::data.table(
  sample_id = rownames(data),
  case_control = syn_data$group
)

### expected data --------------------------------------------------------------

expected_modules <- readRDS("./test_data/coremo_modules.rds")
expected_hvg_data <- readRDS("./test_data/hvg_data.rds")

## run coremo ------------------------------------------------------------------

cor_test <- bulk_coexp(raw_data = data, meta_data = meta_data) %>%
  preprocess_bulk_coexp(., hvg = 0.5, .verbose = FALSE) %>%
  cor_module_processing(., cor_method = "spearman", .verbose = FALSE) %>%
  cor_module_check_epsilon(., rbf_func = "gaussian", .verbose = FALSE) %>%
  cor_module_coremo_clustering(
    .,
    coremo_params = params_coremo(epsilon = 1),
    .verbose = FALSE
  ) %>%
  cor_module_coremo_stability(., .verbose = FALSE)

## tests -----------------------------------------------------------------------

### dimensions -----------------------------------------------------------------

cor_matrix <- S7::prop(cor_test, "processed_data")[[
  "correlation_res"
]]$get_cor_matrix(
  .verbose = FALSE
)

expect_equal(
  current = dim(cor_matrix),
  target = c(500, 500),
  info = paste(
    "CoReMo: cor matrix dimensions"
  )
)

### hvg data -------------------------------------------------------------------

hvg_data <- S7::prop(cor_test, "processed_data")[["feature_meta"]]

expect_equal(
  current = hvg_data,
  target = expected_hvg_data,
  info = paste(
    "CoReMo: HVG data"
  )
)

### modules --------------------------------------------------------------------

module_data <- S7::prop(cor_test, "outputs")[["final_modules"]]

expect_equal(
  current = module_data,
  target = expected_modules,
  info = paste(
    "CoReMo: expected modules"
  )
)
