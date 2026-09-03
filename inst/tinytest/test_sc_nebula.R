# NEBULA tests -----------------------------------------------------------------

source("helper_sc.R", local = TRUE)

library(magrittr)

test_temp_dir <- sc_test_dir("sc_nebula")

## test parameters -------------------------------------------------------------

min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
hvg_to_keep <- 50L
no_pcs <- 20L
n_samples <- 6L
n_genes_test <- 25L

## synthetic test data ---------------------------------------------------------

single_cell_test_data <- sc_test_fixture(
  min_lib_size = min_lib_size,
  min_genes_exp = min_genes_exp,
  min_cells_exp = min_cells_exp,
  hvg_to_keep = hvg_to_keep,
  no_pcs = no_pcs,
  syn_data_params = params_sc_synthetic_data(
    n_samples = n_samples,
    sample_bias = "even"
  )
)

# sample_id is the subject the cells came from, condition a per-subject label
single_cell_test_data$obs[,
  condition := ifelse(
    sample_id %in% c("sample_1", "sample_2", "sample_3"),
    "ctr",
    "trt"
  )
]

## object gen ------------------------------------------------------------------

sc_object <- sc_test_object(
  test_temp_dir,
  single_cell_test_data,
  sc_qc_param = sc_test_qc_params(single_cell_test_data, target_size = 1000)
)

genes_to_test <- head(get_gene_names(sc_object), n_genes_test)

# tests ------------------------------------------------------------------------

## parameters ------------------------------------------------------------------

nebula_params <- params_nebula()

expect_true(
  current = bixverse:::checkNebulaParams(nebula_params),
  info = "the default NEBULA params pass their own check"
)

expect_equal(
  current = nebula_params$nebula_method,
  target = "ln",
  info = "NEBULA params default to the LN method"
)

expect_error(
  current = params_nebula(nebula_method = "nonsense"),
  info = "params_nebula rejects an unknown method"
)

expect_error(
  current = params_nebula(min_sigma = 20, max_sigma = 10),
  pattern = "below",
  info = "params_nebula rejects min_sigma above max_sigma"
)

expect_error(
  current = params_nebula(min_phi = 20, max_phi = 10),
  pattern = "below",
  info = "params_nebula rejects min_phi above max_phi"
)

expect_error(
  current = params_nebula(gene_batch_size = 0L),
  info = "params_nebula rejects a zero gene batch size"
)

expect_true(
  current = is.character(bixverse:::checkNebulaParams(
    utils::modifyList(nebula_params, list(min_sigma = 100))
  )),
  info = "checkNebulaParams catches bounds that cross"
)

expect_true(
  current = is.character(bixverse:::checkNebulaParams(nebula_params[-1])),
  info = "checkNebulaParams returns a message when an element is missing"
)

## input validation ------------------------------------------------------------

expect_error(
  current = nebula_sc(
    object = sc_object,
    subject_col = "not_a_column",
    design = ~condition,
    genes_to_use = genes_to_test,
    .verbose = FALSE
  ),
  info = "nebula_sc rejects a subject column that is not in the obs table"
)

expect_error(
  current = nebula_sc(
    object = sc_object,
    subject_col = "sample_id",
    design = ~not_a_column,
    genes_to_use = genes_to_test,
    .verbose = FALSE
  ),
  info = "nebula_sc rejects a design variable that is not in the obs table"
)

expect_error(
  current = nebula_sc(
    object = sc_object,
    subject_col = "sample_id",
    design = ~condition,
    coef = 2L,
    contrast = c(0, 1),
    genes_to_use = genes_to_test,
    .verbose = FALSE
  ),
  pattern = "not both",
  info = "nebula_sc rejects coef and contrast together"
)

expect_error(
  current = nebula_sc(
    object = sc_object,
    subject_col = "sample_id",
    design = ~condition,
    contrast = c(0, 1, 0),
    genes_to_use = genes_to_test,
    .verbose = FALSE
  ),
  info = "nebula_sc rejects a contrast of the wrong length"
)

expect_error(
  current = nebula_sc(
    object = sc_object,
    subject_col = "sample_id",
    design = ~condition,
    genes_to_use = genes_to_test,
    offset = c(1, 2, 3),
    .verbose = FALSE
  ),
  info = "nebula_sc rejects an offset that does not match the cell count"
)

## single cell results ---------------------------------------------------------

nebula_res <- nebula_sc(
  object = sc_object,
  subject_col = "sample_id",
  design = ~condition,
  genes_to_use = genes_to_test,
  .verbose = FALSE
)

expect_inherits(
  current = nebula_res,
  class = "ScNebula",
  info = "nebula_sc returns the ScNebula class"
)

expect_true(
  current = checkmate::testDataTable(nebula_res$results),
  info = "the NEBULA results are a data.table"
)

expect_true(
  current = checkmate::testNames(
    names(nebula_res$results),
    must.include = c(
      "gene_id",
      "log_fc",
      "effect_se",
      "z",
      "p_value",
      "fdr",
      "subject_overdispersion",
      "cell_overdispersion",
      "convergence",
      "sigma_at_bound"
    )
  ),
  info = "the NEBULA results carry the expected columns"
)

expect_true(
  current = "cell_overdispersion_shrunk" %in% names(nebula_res$results),
  info = "the shrunk overdispersion is present when shrinkage is on"
)

expect_true(
  current = all(nebula_res$results$gene_id %in% genes_to_test),
  info = "only the requested genes come back"
)

expect_equal(
  current = colnames(nebula_res$coefficients),
  target = c("(Intercept)", "conditiontrt"),
  info = "the coefficient matrix carries the design column names"
)

expect_equal(
  current = dim(nebula_res$coefficients),
  target = c(nrow(nebula_res$results), 2L),
  info = "the coefficient matrix is genes x coefficients"
)

expect_equal(
  current = dim(nebula_res$se),
  target = dim(nebula_res$coefficients),
  info = "the standard errors match the coefficients in shape"
)

expect_equal(
  current = rownames(nebula_res$coefficients),
  target = nebula_res$results$gene_id,
  info = "the coefficient rows line up with the results table"
)

## the statistics are internally consistent ------------------------------------

expect_true(
  current = all(
    nebula_res$results$p_value >= 0 & nebula_res$results$p_value <= 1
  ),
  info = "p-values are in [0, 1]"
)

expect_true(
  current = all(nebula_res$results$fdr >= nebula_res$results$p_value - 1e-12),
  info = "the FDR is never below the raw p-value"
)

expect_true(
  current = all(nebula_res$results$effect_se > 0),
  info = "the standard errors are positive"
)

# the Wald statistic is the effect over its standard error, by definition
expect_equal(
  current = nebula_res$results$z,
  target = nebula_res$results$log_fc / nebula_res$results$effect_se,
  tolerance = 1e-8,
  info = "the Wald statistic is log_fc / effect_se"
)

# and the p-value is its two-sided normal tail
expect_equal(
  current = nebula_res$results$p_value,
  target = 2 * stats::pnorm(-abs(nebula_res$results$z)),
  tolerance = 1e-8,
  info = "the p-value is the two-sided normal tail of the Wald statistic"
)

expect_true(
  current = all(nebula_res$results$cell_overdispersion > 0),
  info = "the cell-level overdispersions are positive"
)

expect_true(
  current = all(nebula_res$results$subject_overdispersion >= 0),
  info = "the subject-level overdispersions are non-negative"
)

expect_true(
  current = all(nebula_res$results$convergence > -20L),
  info = "every gene converged on this fixture"
)

## the params are recorded -----------------------------------------------------

params <- get_params(nebula_res)

expect_equal(
  current = params$subject_col,
  target = "sample_id",
  info = "the subject column is recorded in the params"
)

expect_equal(
  current = params$n_subjects,
  target = n_samples,
  info = "the number of subjects is recorded in the params"
)

expect_equal(
  current = params$tested,
  target = "conditiontrt",
  info = "the tested coefficient is recorded by name"
)

expect_true(
  current = checkmate::testString(get_params(nebula_res, to_json = TRUE)),
  info = "the params serialise to JSON"
)

## coefficient selection -------------------------------------------------------

nebula_by_name <- nebula_sc(
  object = sc_object,
  subject_col = "sample_id",
  design = ~condition,
  coef = "conditiontrt",
  genes_to_use = genes_to_test,
  .verbose = FALSE
)

expect_equal(
  current = nebula_by_name$results$log_fc,
  target = nebula_res$results$log_fc,
  info = "naming the coefficient matches the default last-column choice"
)

nebula_contrast <- nebula_sc(
  object = sc_object,
  subject_col = "sample_id",
  design = ~condition,
  contrast = c(0, 1),
  genes_to_use = genes_to_test,
  .verbose = FALSE
)

expect_equal(
  current = nebula_contrast$results$log_fc,
  target = nebula_res$results$log_fc,
  tolerance = 1e-8,
  info = "a unit contrast reproduces testing that coefficient directly"
)

## the fit does not depend on the cell ordering --------------------------------

# NEBULA sorts the cells by subject internally, so a reordered obs table has to
# give the same answer
nebula_shrunk_off <- nebula_sc(
  object = sc_object,
  subject_col = "sample_id",
  design = ~condition,
  genes_to_use = genes_to_test,
  nebula_params = params_nebula(shrink_dispersion = FALSE),
  .verbose = FALSE
)

expect_false(
  current = "cell_overdispersion_shrunk" %in%
    names(nebula_shrunk_off$results),
  info = "no shrunk overdispersion column when shrinkage is off"
)

expect_equal(
  current = nebula_shrunk_off$results$log_fc,
  target = nebula_res$results$log_fc,
  info = "the shrinkage does not move the fixed effects"
)

## the expression filter bites -------------------------------------------------

nebula_strict <- nebula_sc(
  object = sc_object,
  subject_col = "sample_id",
  design = ~condition,
  genes_to_use = genes_to_test,
  nebula_params = params_nebula(cpc = 5),
  .verbose = FALSE
)

expect_true(
  current = nrow(nebula_strict$results) < nrow(nebula_res$results),
  info = "a higher counts-per-cell cutoff drops genes"
)

# losing every gene is an error rather than an empty table, so the caller
# cannot mistake a wiped-out filter for a null result
expect_error(
  current = nebula_sc(
    object = sc_object,
    subject_col = "sample_id",
    design = ~condition,
    genes_to_use = genes_to_test,
    nebula_params = params_nebula(cpc = 1e6),
    .verbose = FALSE
  ),
  pattern = "expression filter",
  info = "NEBULA errors when the filter leaves nothing"
)

## meta cells ------------------------------------------------------------------

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)
sc_object <- calculate_pca_sc(sc_object, no_pcs = no_pcs, .verbose = FALSE)
sc_object <- find_neighbours_sc(sc_object, .verbose = FALSE)

mc_object <- generate_bt_meta_cells_sc(
  sc_object,
  sc_meta_cell_params = params_sc_bt_metacells(target_no_metacells = 120L),
  .verbose = FALSE
)

# the meta cell obs table carries no donor labels, so map them over from the
# originating cells by majority vote
sc_obs <- sc_object[[c("cell_idx", "sample_id")]]
mc_object[["sample_id"]] <- vapply(
  get_sc_obs(mc_object)$original_cell_idx,
  \(idx) {
    tab <- table(sc_obs$sample_id[match(idx, sc_obs$cell_idx)])
    names(tab)[which.max(tab)]
  },
  character(1)
)
mc_object[["condition"]] <- ifelse(
  get_sc_obs(mc_object)$sample_id %in%
    c("sample_1", "sample_2", "sample_3"),
  "ctr",
  "trt"
)

mc_genes <- head(S7::prop(mc_object, "var_table")$gene_id, n_genes_test)

nebula_mc_res <- nebula_mc(
  object = mc_object,
  subject_col = "sample_id",
  design = ~condition,
  genes_to_use = mc_genes,
  .verbose = FALSE
)

expect_inherits(
  current = nebula_mc_res,
  class = "ScNebula",
  info = "nebula_mc returns the ScNebula class"
)

expect_true(
  current = checkmate::testNames(
    names(nebula_mc_res$results),
    must.include = c("gene_id", "log_fc", "z", "p_value", "fdr")
  ),
  info = "the meta cell NEBULA results carry the expected columns"
)

expect_equal(
  current = colnames(nebula_mc_res$coefficients),
  target = c("(Intercept)", "conditiontrt"),
  info = "the meta cell coefficient matrix carries the design column names"
)

expect_true(
  current = all(
    nebula_mc_res$results$p_value >= 0 & nebula_mc_res$results$p_value <= 1
  ),
  info = "meta cell p-values are in [0, 1]"
)

expect_equal(
  current = nebula_mc_res$results$z,
  target = nebula_mc_res$results$log_fc / nebula_mc_res$results$effect_se,
  tolerance = 1e-8,
  info = "the meta cell Wald statistic is log_fc / effect_se"
)

expect_equal(
  current = get_params(nebula_mc_res)$n_subjects,
  target = n_samples,
  info = "the meta cell run sees the same number of subjects"
)

# aggregation smooths within a subject, so the cell-level term shrinks. That is
# the whole reason the two runs are not comparable.
expect_true(
  current = stats::median(nebula_mc_res$results$cell_overdispersion) <
    stats::median(nebula_res$results$cell_overdispersion),
  info = "the meta cell overdispersion is below the single cell one"
)

expect_error(
  current = nebula_mc(
    object = mc_object,
    subject_col = "not_a_column",
    design = ~condition,
    genes_to_use = mc_genes,
    .verbose = FALSE
  ),
  info = "nebula_mc rejects a subject column that is not in the obs table"
)


## subsets ---------------------------------------------------------------------

# Running NEBULA within a cell type is the normal usage, so the subset path
# matters. A subset shares its parent's counts file and carries only its own
# obs rows, with cell_idx still in the parent's index space, so the subject
# vector has to be sized off the store rather than off the obs table.
sc_object[["cell_grp_chr"]] <- as.character(sc_object[[]]$cell_grp)
first_grp <- sort(unique(sc_object[[]]$cell_grp_chr))[1]

sc_subset <- SingleCellsSubset(
  sc_object,
  grouping_column = "cell_grp_chr",
  group = first_grp
)

expect_true(
  current = S7::prop(sc_subset, "dims")[1] < S7::prop(sc_object, "dims")[1],
  info = "the subset holds fewer cells than the parent"
)

nebula_subset <- nebula_sc(
  object = sc_subset,
  subject_col = "sample_id",
  design = ~condition,
  genes_to_use = genes_to_test,
  .verbose = FALSE
)

expect_inherits(
  current = nebula_subset,
  class = "ScNebula",
  info = "nebula_sc runs on a SingleCellsSubset"
)

expect_equal(
  current = get_params(nebula_subset)$n_cells,
  target = S7::prop(sc_subset, "dims")[1],
  info = "the subset run sees exactly the subset's cells"
)

expect_true(
  current = all(
    nebula_subset$results$p_value >= 0 & nebula_subset$results$p_value <= 1
  ),
  info = "the subset run produces usable p-values"
)
# clean up ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir)
