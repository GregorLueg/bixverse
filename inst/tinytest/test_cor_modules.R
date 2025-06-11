# cor module tests -------------------------------------------------------------

library(magrittr)

## data ------------------------------------------------------------------------

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
