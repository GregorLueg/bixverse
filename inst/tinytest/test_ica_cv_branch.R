# ICA cross-validate branch lock-in --------------------------------------------

library(magrittr)

## data ------------------------------------------------------------------------

synth <- synthetic_bulk_cor_matrix(
  num_samples = 30L,
  num_genes = 150L,
  add_modules = TRUE,
  module_sizes = c(30L, 30L),
  seed = 42L
)

counts <- synth$counts
meta_data <- data.table::data.table(sample_id = colnames(counts))

## pipeline with cross_validate = TRUE -----------------------------------------

ica_cv_test <- BulkCoExp(
  raw_data = t(counts),
  meta_data = meta_data
) %>%
  preprocess_bulk_coexp(
    scaling = TRUE,
    scaling_type = "normal",
    hvg = NULL,
    .verbose = FALSE
  ) %>%
  ica_processing(.verbose = FALSE)

expect_silent_run <- try(
  ica_cv_test <- ica_stabilised_results(
    ica_cv_test,
    no_comp = 3L,
    ica_type = "logcosh",
    iter_params = params_ica_randomisation(
      cross_validate = TRUE,
      random_init = 3L,
      folds = 3L
    ),
    random_seed = 42L,
    .verbose = FALSE
  ),
  silent = TRUE
)

expect_false(
  current = inherits(expect_silent_run, "try-error"),
  info = "ICA CV branch runs without error"
)

result <- get_results(ica_cv_test)

expect_true(
  current = inherits(result, "BulkModuleResult"),
  info = "ICA CV branch returns BulkModuleResult"
)

expect_true(
  current = data.table::is.data.table(get_modules(result)),
  info = "ICA CV branch - get_modules() returns a data.table"
)
