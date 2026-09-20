# residual-based single cell path ----------------------------------------------

# scTransform and analytic Pearson, the residual HVG selection they feed and the
# residual PCA on top. The numerical parity against the R sctransform package is
# pinned in bixverse-rs; what is checked here is the marshalling, the guards and
# the object state.

library(bixverse)

source("helper_sc.R", local = TRUE)

set.seed(123L)

## setup -----------------------------------------------------------------------

fixture <- sc_test_fixture()
test_dir <- sc_test_dir("sc_residuals")
sc_object <- sc_test_object(dir = test_dir, fixture = fixture)

n_cells <- length(get_cells_to_keep(sc_object))
n_genes <- dim(sc_object)[2]

# the synthetic obs carries three cell types, which stand in for samples
obs <- get_sc_obs(sc_object, filtered = TRUE)
group_col <- "cell_grp"

## analytic Pearson ------------------------------------------------------------

### fitting --------------------------------------------------------------------

sc_apr <- fit_residuals_sc(
  object = sc_object,
  method = "analytic_pearson",
  .verbose = FALSE
)

apr_fit <- get_residual_fit(sc_apr)

expect_inherits(
  current = apr_fit,
  class = "ScResidualFit",
  info = "The fit comes back as its own S3 class"
)

expect_equal(
  current = apr_fit$method,
  target = "analytic_pearson",
  info = "The fit records the method it was built with"
)

expect_equal(
  current = apr_fit$n_groups,
  target = 1L,
  info = "An ungrouped fit has a single model"
)

expect_equal(
  current = length(apr_fit$cell_indices),
  target = n_cells,
  info = "The fit covers every kept cell"
)

expect_true(
  current = !is.unsorted(apr_fit$genes),
  info = "The gene axis is ascending, which the binary searches rely on"
)

expect_true(
  current = all(apr_fit$genes >= 0) && all(apr_fit$genes < n_genes),
  info = "The gene axis stays inside the store"
)

expect_true(
  current = all(apr_fit$cell_totals > 0),
  info = "Per-cell totals over the retained genes are positive"
)

expect_true(
  current = is.null(apr_fit$log10_umi),
  info = "The analytic Pearson model carries no scTransform offset"
)

### residuals against a hand computation ---------------------------------------

# the analytic Pearson residual is closed form, so R can produce the reference
# directly:  mu = gene_sum * cell_total / total,
#            r  = (y - mu) / sqrt(mu + mu^2 / theta),  clipped at +/- sqrt(n)
raw_counts <- sc_object[,
  as.integer(apr_fit$genes) + 1L,
  assay = "raw",
  return_format = "gene",
  use_cells_to_keep = TRUE
]
raw_counts <- as.matrix(raw_counts)

apr_model <- apr_fit$models[[1]]
mu_r <- outer(apr_fit$cell_totals, apr_model$gene_sums) / apr_model$total
resid_r <- (raw_counts - mu_r) / sqrt(mu_r + mu_r^2 / apr_model$theta)
resid_r <- pmin(pmax(resid_r, apr_model$clip_min), apr_model$clip_max)

# sample variance, matching the (n - 1) denominator the Rust reducer uses
var_r <- apply(resid_r, 2, var)

sc_apr <- find_hvg_sc(
  object = sc_apr,
  hvg_no = fixture$hvg_to_keep,
  hvg_params = params_sc_hvg(method = "residual"),
  .verbose = FALSE
)

var_rust <- get_sc_var(sc_apr)$residual_variance[
  as.integer(apr_fit$genes) + 1L
]

expect_equivalent(
  current = var_rust,
  target = var_r,
  tolerance = 1e-5,
  info = "Residual variance matches the closed form computed in R"
)

### hvg selection --------------------------------------------------------------

expect_equal(
  current = length(get_hvg(sc_apr)),
  target = fixture$hvg_to_keep,
  info = "A single group selects exactly hvg_no features"
)

expect_true(
  current = all(get_hvg(sc_apr) %in% apr_fit$genes),
  info = "Every selected feature is one the model covers"
)

hvg_r <- order(var_r, decreasing = TRUE)[1:fixture$hvg_to_keep]
expect_equal(
  current = sort(get_hvg(sc_apr)),
  target = sort(as.integer(apr_fit$genes[hvg_r])),
  info = "The selection is the top N by residual variance"
)

### variance column ------------------------------------------------------------

var_table <- get_sc_var(sc_apr)

expect_true(
  current = "residual_variance" %in% names(var_table),
  info = "The residual variance lands in the var table"
)

expect_equal(
  current = nrow(var_table),
  target = n_genes,
  info = "The var table keeps one row per gene in the store"
)

expect_equal(
  current = sum(!is.na(var_table$residual_variance)),
  target = length(apr_fit$genes),
  info = "Genes outside the model are NA, not zero"
)

## scTransform -----------------------------------------------------------------

sc_sct <- fit_residuals_sc(
  object = sc_object,
  method = "sctransform",
  .verbose = FALSE
)

sct_fit <- get_residual_fit(sc_sct)

expect_equal(
  current = sct_fit$method,
  target = "sctransform",
  info = "The scTransform fit records its method"
)

sct_model <- sct_fit$models[[1]]

expect_equal(
  current = sct_model$n_coef,
  target = 1L,
  info = "Without covariates the design is the intercept alone"
)

expect_equal(
  current = length(sct_model$theta),
  target = length(sct_fit$genes),
  info = "One dispersion per modelled gene"
)

expect_equal(
  current = length(sct_model$coefficients),
  target = length(sct_fit$genes) * sct_model$n_coef,
  info = "The coefficient block is genes by coefficients"
)

expect_equal(
  current = length(sct_fit$log10_umi),
  target = n_cells,
  info = "The offset carries one entry per selected cell"
)

expect_true(
  current = all(is.finite(sct_model$theta) | is.infinite(sct_model$theta)),
  info = "Dispersions are finite or the Poisson limit, never NaN"
)

sc_sct <- find_hvg_sc(
  object = sc_sct,
  hvg_no = fixture$hvg_to_keep,
  hvg_params = params_sc_hvg(method = "residual"),
  .verbose = FALSE
)

# what the transform is for: a gene with no structure beyond its mean should
# come out with a residual variance near one
sct_variance <- get_sc_var(sc_sct)$residual_variance

expect_true(
  current = abs(median(sct_variance, na.rm = TRUE) - 1) < 0.5,
  info = "Median residual variance sits near one"
)

expect_true(
  current = all(sct_variance[!is.na(sct_variance)] >= 0),
  info = "A variance is never negative"
)

### covariates -----------------------------------------------------------------

sc_cov <- fit_residuals_sc(
  object = sc_object,
  method = "sctransform",
  covariate_columns = "nnz",
  .verbose = FALSE
)

cov_fit <- get_residual_fit(sc_cov)

expect_equal(
  current = cov_fit$models[[1]]$n_coef,
  target = 2L,
  info = "One covariate adds a column to the intercept-only design"
)

expect_equal(
  current = cov_fit$covariate_names,
  target = "nnz",
  info = "The covariate order is recorded so it can be checked on reuse"
)

expect_error(
  current = fit_residuals_sc(
    object = sc_object,
    method = "analytic_pearson",
    covariate_columns = "nnz",
    .verbose = FALSE
  ),
  pattern = "only supported for method",
  info = "The analytic Pearson model has no design to add covariates to"
)

expect_error(
  current = fit_residuals_sc(
    object = sc_object,
    method = "sctransform",
    covariate_columns = "not_a_column",
    .verbose = FALSE
  ),
  pattern = "not found in the observation table",
  info = "A covariate column that does not exist is named in the error"
)

expect_error(
  current = fit_residuals_sc(
    object = sc_object,
    method = "sctransform",
    covariate_columns = group_col,
    .verbose = FALSE
  ),
  pattern = "not numeric",
  info = "A factor covariate is refused rather than dummy coded silently"
)

## grouped fits ----------------------------------------------------------------

sc_grouped <- fit_residuals_sc(
  object = sc_object,
  method = "analytic_pearson",
  group_column = group_col,
  .verbose = FALSE
)

grouped_fit <- get_residual_fit(sc_grouped)
n_groups <- length(unique(obs[[group_col]]))

expect_equal(
  current = grouped_fit$n_groups,
  target = n_groups,
  info = "One model per level of the grouping column"
)

expect_equal(
  current = length(grouped_fit$models),
  target = n_groups,
  info = "The model list matches the group count"
)

expect_equal(
  current = sort(unique(grouped_fit$group_of_cell)),
  target = seq_len(n_groups) - 1L,
  info = "Group labels densely cover 0 to n_groups - 1"
)

expect_true(
  current = length(grouped_fit$genes) <=
    min(purrr::map_int(grouped_fit$models, \(m) length(m$genes))),
  info = "The shared axis is the intersection of the per-group sets"
)

expect_equal(
  current = grouped_fit$group_column,
  target = group_col,
  info = "The fit remembers what it was grouped on"
)

### the union rule -------------------------------------------------------------

hvg_per_group <- 10L
sc_grouped <- find_hvg_sc(
  object = sc_grouped,
  hvg_no = hvg_per_group,
  hvg_params = params_sc_hvg(method = "residual"),
  .verbose = FALSE
)

n_selected <- length(get_hvg(sc_grouped))

expect_true(
  current = n_selected >= hvg_per_group &&
    n_selected <= hvg_per_group * n_groups,
  info = paste(
    "hvg_no is a per-group count and the groups are unioned, so the result",
    "sits between hvg_no and n_groups * hvg_no"
  )
)

## residual pca ----------------------------------------------------------------

### the guards -----------------------------------------------------------------

expect_error(
  current = calculate_pca_sc(
    object = sc_apr,
    no_pcs = fixture$no_pcs,
    residuals = TRUE,
    .verbose = FALSE
  ),
  pattern = "normalise_variance",
  info = "Variance normalisation is refused, not silently overridden"
)

expect_error(
  current = calculate_pca_sc(
    object = sc_apr,
    no_pcs = fixture$no_pcs,
    pca_params = params_sc_pca(normalise_variance = FALSE, clr = TRUE),
    residuals = TRUE,
    .verbose = FALSE
  ),
  pattern = "PFlogPF",
  info = "The CLR transform belongs to the normalised layer"
)

expect_error(
  current = calculate_pca_sc(
    object = sc_apr,
    no_pcs = fixture$no_pcs,
    pca_params = params_sc_pca(normalise_variance = FALSE),
    sparse_svd = TRUE,
    residuals = TRUE,
    .verbose = FALSE
  ),
  pattern = "no sparse solver",
  info = "A residual column is dense, so there is no sparse path"
)

expect_error(
  current = calculate_pca_sc(
    object = sc_object,
    no_pcs = fixture$no_pcs,
    pca_params = params_sc_pca(normalise_variance = FALSE),
    hvg = fixture$genes_pass[1:10],
    residuals = TRUE,
    .verbose = FALSE
  ),
  pattern = "No fitted residual model",
  info = "The PCA asserts on a missing fit rather than warning"
)

### the happy path -------------------------------------------------------------

sc_apr <- calculate_pca_sc(
  object = sc_apr,
  no_pcs = fixture$no_pcs,
  pca_params = params_sc_pca(normalise_variance = FALSE),
  residuals = TRUE,
  .verbose = FALSE
)

pca_factors <- get_pca_factors(sc_apr)

expect_equal(
  current = dim(pca_factors),
  target = c(n_cells, fixture$no_pcs),
  info = "The scores are cells by components"
)

expect_true(
  current = !any(is.nan(pca_factors)) && !any(is.na(pca_factors)),
  info = "No NaN leaks out of the residual arithmetic"
)

expect_equal(
  current = nrow(get_pca_loadings(sc_apr)),
  target = length(get_hvg(sc_apr)),
  info = "The loadings cover the selected features"
)

expect_equal(
  current = length(get_pca_singular_val(sc_apr)),
  target = fixture$no_pcs,
  info = "One singular value per requested component"
)

## staleness -------------------------------------------------------------------

sc_stale <- suppressWarnings(set_cells_to_keep(
  sc_apr,
  get_cells_to_keep(sc_apr)[1:100] + 1L
))

expect_error(
  current = find_hvg_sc(
    object = sc_stale,
    hvg_no = 10L,
    hvg_params = params_sc_hvg(method = "residual"),
    .verbose = FALSE
  ),
  pattern = "cell filter moved",
  info = "A fit is refused once the cell selection it was built on changed"
)

## unsupported combinations ----------------------------------------------------

expect_error(
  current = get_hvg_data_sc(
    object = sc_apr,
    hvg_params = params_sc_hvg(method = "residual"),
    .verbose = FALSE
  ),
  pattern = "does not support method",
  info = "The non-mutating getter would have to refit, so it refuses"
)

expect_error(
  current = find_hvg_batch_aware_sc(
    object = sc_apr,
    batch_column = group_col,
    hvg_params = params_sc_hvg(method = "residual"),
    .verbose = FALSE
  ),
  pattern = "does not support method",
  info = "Per-group residual selection happens at fitting time instead"
)

expect_error(
  current = sct_corrected_counts_sc(sc_apr, .verbose = FALSE),
  pattern = "need a scTransform fit",
  info = "The analytic Pearson model has no corrected-count equivalent"
)

## parameter validation --------------------------------------------------------

expect_error(
  current = params_sc_sctransform(clip_min = -5),
  pattern = "must be given together",
  info = "Half a clipping range would silently become no range at all"
)

expect_error(
  current = params_sc_apr(clip_min = 5, clip_max = -5),
  pattern = "smaller than",
  info = "An inverted clipping range is caught in R"
)

expect_silent(
  current = params_sc_apr(theta = Inf),
  info = "Infinite theta is the Poisson limit and a legitimate choice"
)

expect_error(
  current = params_sc_apr(theta = 0),
  info = "A non-positive dispersion is not"
)

expect_silent(
  current = params_sc_apr(min_cells = 0L),
  info = "Zero min_cells keeps every gene, which is allowed"
)

## meta cells ------------------------------------------------------------------

# meta cell counts are summed UMIs, so the negative binomial still applies. The
# path runs in memory through the same Rust code behind a different reader.

mc_object <- generate_bt_meta_cells_sc(
  object = sc_test_prepped(sc_object, fixture),
  sc_meta_cell_params = params_sc_bt_metacells(target_no_metacells = 50L),
  .verbose = FALSE
)

mc_apr <- fit_residuals_sc(
  object = mc_object,
  method = "analytic_pearson",
  residual_params = params_sc_apr(min_cells = 2L),
  .verbose = FALSE
)

mc_fit <- get_residual_fit(mc_apr)

expect_inherits(
  current = mc_fit,
  class = "ScResidualFit",
  info = "Meta cells carry the same fit class as the on-disk classes"
)

expect_equal(
  current = length(mc_fit$cell_indices),
  target = nrow(S7::prop(mc_object, "obs_table")),
  info = "The fit covers every meta cell"
)

mc_apr <- find_hvg_sc(
  object = mc_apr,
  hvg_no = 15L,
  hvg_params = params_sc_hvg(method = "residual"),
  .verbose = FALSE
)

mc_hvg <- get_hvg(mc_apr)

expect_equal(
  current = length(mc_hvg),
  target = 15L,
  info = "A single group selects exactly hvg_no features"
)

expect_true(
  current = min(mc_hvg) >= 1L,
  info = "Meta cells index genes 1-based, unlike the on-disk classes"
)

mc_apr <- calculate_pca_sc(
  object = mc_apr,
  no_pcs = 5L,
  pca_params = params_sc_pca(normalise_variance = FALSE),
  residuals = TRUE,
  .verbose = FALSE
)

expect_equal(
  current = dim(get_pca_factors(mc_apr)),
  target = c(nrow(S7::prop(mc_object, "obs_table")), 5L),
  info = "The meta cell scores are meta cells by components"
)

expect_true(
  current = !any(is.nan(get_pca_factors(mc_apr))),
  info = "No NaN leaks out of the in-memory residual arithmetic"
)

expect_error(
  current = calculate_pca_sc(
    object = mc_apr,
    no_pcs = 5L,
    residuals = TRUE,
    .verbose = FALSE
  ),
  pattern = "normalise_variance",
  info = "The meta cell path enforces the same PCA settings"
)

expect_error(
  current = find_hvg_sc(
    object = mc_object,
    hvg_no = 5L,
    hvg_params = params_sc_hvg(method = "residual"),
    .verbose = FALSE
  ),
  pattern = "No fitted residual model",
  info = "Meta cells assert on a missing fit too"
)

## persistence and shape migration ---------------------------------------------

save_sc_exp_to_disk(sc_apr, type = "qs2")

reloaded <- load_existing(SingleCells(dir_data = test_dir), .verbose = FALSE)
reloaded_fit <- get_residual_fit(reloaded)

expect_equal(
  current = reloaded_fit$method,
  target = apr_fit$method,
  info = "The fitted model survives a save and load round trip"
)

expect_equal(
  current = reloaded_fit$genes,
  target = apr_fit$genes,
  info = "The gene axis comes back unchanged"
)

expect_silent(
  current = find_hvg_sc(
    object = reloaded,
    hvg_no = fixture$hvg_to_keep,
    hvg_params = params_sc_hvg(method = "residual"),
    .verbose = FALSE
  ),
  info = "A reloaded fit is still usable, so its stamp survived too"
)

### objects saved before the slot existed --------------------------------------

memory_path <- file.path(test_dir, "memory.qs2")
saved <- qs2::qs_read(memory_path)

# what a file written by an older bixverse looks like
saved$sc_cache[["residual_fit"]] <- NULL
qs2::qs_save(saved, memory_path)

old_object <- load_existing(SingleCells(dir_data = test_dir), .verbose = FALSE)
old_cache <- get_sc_cache(old_object)

expect_equal(
  current = sort(names(old_cache)),
  target = sort(names(bixverse:::new_sc_cache())),
  info = "A cache saved without the slot is brought up to the current shape"
)

expect_true(
  current = is.null(old_cache$residual_fit),
  info = "The missing slot comes back empty rather than absent"
)

expect_false(
  current = is.null(old_cache$pca_factors),
  info = "Migration keeps everything the saved cache did carry"
)

### slots the constructor no longer knows --------------------------------------

saved$sc_cache[["legacy_slot"]] <- "dropped"
qs2::qs_save(saved, memory_path)

expect_warning(
  current = load_existing(SingleCells(dir_data = test_dir), .verbose = FALSE),
  pattern = "legacy_slot",
  info = "An unknown slot is named on the way out, not silently kept"
)

## cleanup ---------------------------------------------------------------------

sc_test_cleanup(test_dir)
