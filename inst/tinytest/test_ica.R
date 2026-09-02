# tests ica --------------------------------------------------------------------

S <- cbind(sin((1:1000) / 20), rep((((1:200) - 100) / 100), 5))
A <- matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
X <- S %*% A

## ica -------------------------------------------------------------------------

### whitening ------------------------------------------------------------------

expected_k <- matrix(c(-0.8243566, -2.2245938, -1.7226241, 1.0645727), 2, 2)

# prepare the whitening
whitening_res <- rs_prepare_whitening(
  x = X,
  fast_svd = TRUE,
  seed = 42L,
  rank = NULL,
  oversampling = NULL,
  n_power_iter = NULL
)

expect_equal(
  current = whitening_res$k,
  target = expected_k,
  tolerance = 10e-6,
  info = paste(
    "ICA pre-whitening matrix k"
  )
)

### logcosh implementation -----------------------------------------------------

# conscious complete randomisation
random_seed <- as.integer(sample(1:1000, 1))

# logcosh implementation
rs_ica_res_logcosh <- fast_ica_rust_helper(
  X = whitening_res$x,
  K = whitening_res$k,
  n_icas = 2L,
  ica_fun = "logcosh",
  seed = random_seed
)

# pending on seed the sign and the position of the extracted signal
# changes. However, in one case the correlation should be close to abs(1),
# indicating that ICA correctly identified the two sources of the signal

logcosh_correlation_signal_11 <- abs(cor(S[, 1], rs_ica_res_logcosh$S[1, ]))
logcosh_correlation_signal_22 <- abs(cor(S[, 2], rs_ica_res_logcosh$S[2, ]))
logcosh_correlation_signal_12 <- abs(cor(S[, 1], rs_ica_res_logcosh$S[2, ]))
logcosh_correlation_signal_21 <- abs(cor(S[, 2], rs_ica_res_logcosh$S[1, ]))

logcosh_correctly_reconstructed <- purrr::map_lgl(
  c(
    logcosh_correlation_signal_11,
    logcosh_correlation_signal_22,
    logcosh_correlation_signal_12,
    logcosh_correlation_signal_21
  ),
  ~ {
    .x > 0.99
  }
)

expect_true(
  sum(logcosh_correctly_reconstructed) == 2,
  info = paste(
    "ICA signal identification with logcosh"
  )
)

#### exp -----------------------------------------------------------------------

# exp implementation
rs_ica_res_exp <- fast_ica_rust_helper(
  X = whitening_res$x,
  K = whitening_res$k,
  n_icas = 2L,
  ica_fun = "exp",
  seed = random_seed
)

exp_correlation_signal_11 <- abs(cor(S[, 1], rs_ica_res_exp$S[1, ]))
exp_correlation_signal_22 <- abs(cor(S[, 2], rs_ica_res_exp$S[2, ]))
exp_correlation_signal_12 <- abs(cor(S[, 1], rs_ica_res_exp$S[2, ]))
exp_correlation_signal_21 <- abs(cor(S[, 2], rs_ica_res_exp$S[1, ]))

exp_correctly_reconstructed <- purrr::map_lgl(
  c(
    logcosh_correlation_signal_11,
    logcosh_correlation_signal_22,
    logcosh_correlation_signal_12,
    logcosh_correlation_signal_21
  ),
  ~ {
    .x > 0.99
  }
)

expect_true(
  sum(exp_correctly_reconstructed) == 2,
  info = paste(
    "ICA signal identification with exp"
  )
)

### full version ---------------------------------------------------------------

#### logcosh -------------------------------------------------------------------

rs_ica_res_logcosh <- fast_ica_rust(
  X = X,
  n_icas = 2L,
  ica_fun = "logcosh",
  seed = random_seed,
  fast_svd = TRUE
)

logcosh_correlation_signal_11 <- abs(cor(S[, 1], rs_ica_res_logcosh$S[1, ]))
logcosh_correlation_signal_22 <- abs(cor(S[, 2], rs_ica_res_logcosh$S[2, ]))
logcosh_correlation_signal_12 <- abs(cor(S[, 1], rs_ica_res_logcosh$S[2, ]))
logcosh_correlation_signal_21 <- abs(cor(S[, 2], rs_ica_res_logcosh$S[1, ]))

logcosh_correctly_reconstructed <- purrr::map_lgl(
  c(
    logcosh_correlation_signal_11,
    logcosh_correlation_signal_22,
    logcosh_correlation_signal_12,
    logcosh_correlation_signal_21
  ),
  ~ {
    .x > 0.99
  }
)

expect_true(
  sum(logcosh_correctly_reconstructed) == 2,
  info = paste(
    "ICA signal identification with logcosh (full)"
  )
)

#### exp -----------------------------------------------------------------------

rs_ica_res_exp <- fast_ica_rust(
  X = X,
  n_icas = 2L,
  ica_fun = "exp",
  seed = random_seed
)

exp_correlation_signal_11 <- abs(cor(S[, 1], rs_ica_res_exp$S[1, ]))
exp_correlation_signal_22 <- abs(cor(S[, 2], rs_ica_res_exp$S[2, ]))
exp_correlation_signal_12 <- abs(cor(S[, 1], rs_ica_res_exp$S[2, ]))
exp_correlation_signal_21 <- abs(cor(S[, 2], rs_ica_res_exp$S[1, ]))

exp_correctly_reconstructed <- purrr::map_lgl(
  c(
    exp_correlation_signal_11,
    exp_correlation_signal_22,
    exp_correlation_signal_12,
    exp_correlation_signal_21
  ),
  ~ {
    .x > 0.99
  }
)

expect_true(
  sum(exp_correctly_reconstructed) == 2,
  info = paste(
    "ICA signal identification with exp (full)"
  )
)

## ica class -------------------------------------------------------------------

### expected data --------------------------------------------------------------

expected_ica_x1_mat <- qs2::qs_read("./test_data/ica_class_x1.qs")
expected_ica_stability_res <- qs2::qs_read(
  "./test_data/ica_class_stability_res.qs"
)
expected_ica_s_logcosh_mat <- qs2::qs_read(
  "./test_data/ica_s_matrix_logcosh.qs"
)
expected_ica_a_logcosh_mat <- qs2::qs_read(
  "./test_data/ica_a_matrix_logcosh.qs"
)
expected_ica_s_exp_mat <- qs2::qs_read(
  "./test_data/ica_s_matrix_exp.qs"
)
expected_ica_a_exp_mat <- qs2::qs_read(
  "./test_data/ica_a_matrix_exp.qs"
)

### gex synthetic data ---------------------------------------------------------

# The non-Gaussian generator is the one ICA is identifiable on, see
# params_synthetic_bulk_rnaseq().
test_data <- synthetic_bulk_cor_matrix(
  params_synthetic_bulk_rnaseq(
    num_samples = 100L,
    num_genes = 1000L,
    module_sizes = c(100L, 100L, 100L),
    generator = "non_gaussian_factor",
    seed = 123L
  )
)

norm_counts <- edgeR::cpm(test_data$counts, log = TRUE)

data <- t(norm_counts)

meta_data <- data.table::data.table(
  sample_id = rownames(data)
)

### run the class (pre processing) ---------------------------------------------

# this stuff is extensively tested in test_cor_modules.R
ica_test <- BulkCoExp(raw_data = data, meta_data = meta_data)

expect_warning(
  current = ica_processing(ica_test, .verbose = TRUE),
  info = "ica bulk coexp - warning without pre-processed"
)

ica_test <- preprocess_bulk_coexp(ica_test, hvg = 0.5, .verbose = FALSE)
ica_test <- ica_processing(ica_test, .verbose = FALSE)

# K is built from the singular vectors of the covariance, and singular vectors
# belonging to near-degenerate values are only defined up to a rotation of
# their subspace. Adjacent singular values here come within 0.04% of each
# other, so an entrywise comparison pins down an arbitrary choice of basis and
# disagrees between BLAS implementations. Assert what K is for instead: it
# whitens, K V K' = I, which no rotation or sign flip can break.
k_mat <- ica_test@processed_data$K
x1_mat <- ica_test@processed_data$X1
cov_mat <- (x1_mat %*% t(x1_mat)) / ncol(x1_mat)

expect_equal(
  current = nrow(k_mat),
  target = nrow(data) - 1L,
  info = "ica bulk coexp - k matrix capped at the achievable rank"
)

expect_true(
  current = max(abs(k_mat %*% cov_mat %*% t(k_mat) - diag(nrow(k_mat)))) < 1e-8,
  info = "ica bulk coexp - k matrix whitens the covariance"
)

# the whitening scale is 1/sqrt(sigma) over descending sigma, so the row norms
# only ever grow. A null-space row would spike instead.
expect_false(
  current = is.unsorted(sqrt(rowSums(k_mat^2))),
  info = "ica bulk coexp - k matrix row norms grow with the component index"
)

expect_equal(
  current = ica_test@processed_data$X1,
  target = expected_ica_x1_mat,
  info = "ica bulk coexp - x1 matrix"
)

if (identical(Sys.getenv("REGEN"), "1")) {
  qs2::qs_save(ica_test@processed_data$X1, "./test_data/ica_class_x1.qs")
}

### ica runs -------------------------------------------------------------------

ica_test <- ica_evaluate_comp(
  ica_test,
  ica_type = "logcosh",
  ncomp_params = params_ica_ncomp(steps = 3L, max_no_comp = 30L),
  .verbose = FALSE
)

# this one should throw a warning, as the loess function looks ugly
expect_warning(
  current = ica_optimal_ncomp(ica_test, .verbose = FALSE, show_plot = FALSE),
  info = paste("ica bulk coexp - warnings from loess")
)

p_individual <- plot_ica_stability_individual(ica_test)
p_ncomp_params <- plot_ica_ncomp_params(ica_test)

expect_true(
  "ggplot" %in% class(p_individual),
  info = paste("ica bulk coexp - individual stability plots")
)

expect_true(
  "ggplot" %in% class(p_ncomp_params),
  info = paste("ica bulk coexp - all parameter plots")
)

ica_stability_res <- get_ica_stability_res(ica_test)

expect_true(
  current = checkmate::test_data_table(ica_stability_res),
  info = paste("ica bulk coexp - stability results class")
)

# seed issues and float precision again
corr_stability <- cor(
  ica_stability_res$median_stability,
  expected_ica_stability_res$median_stability,
  method = "spearman"
)

corr_convergence <- cor(
  ica_stability_res$converged,
  expected_ica_stability_res$converged,
  method = "spearman"
)

corr_mutual_info <- cor(
  ica_stability_res$norm_mutual_information,
  expected_ica_stability_res$norm_mutual_information,
  method = "spearman"
)

expect_true(
  current = corr_stability >= 0.75,
  info = paste("ica bulk coexp - stability: median results")
)

expect_true(
  current = corr_convergence >= 0.95,
  info = paste("ica bulk coexp - stability: covergence results")
)

expect_true(
  current = corr_mutual_info >= 0.95,
  info = paste("ica bulk coexp - stability: mutual info results")
)

if (identical(Sys.getenv("REGEN"), "1")) {
  qs2::qs_save(ica_stability_res, "./test_data/ica_class_stability_res.qs")
}

### class (logcosh) ------------------------------------------------------------

ica_test <- ica_stabilised_results(ica_test, no_comp = 3L, ica_type = "logcosh")

results <- get_results(ica_test)

expect_true(
  current = inherits(results, "BulkModuleResult"),
  info = "ICA class - final_results is a BulkModuleResult"
)

expect_equal(
  current = get_factors(results, which = "gene_loadings"),
  target = t(expected_ica_s_logcosh_mat),
  info = paste("ICA class - gene_loadings (logcosh, transposed S)")
)

expect_equal(
  current = get_factors(results, which = "sample_activity"),
  target = expected_ica_a_logcosh_mat,
  info = paste("ICA class - sample_activity (logcosh, A)")
)

expect_true(
  current = data.table::is.data.table(get_modules(results)),
  info = "ICA class - get_modules() returns a data.table (logcosh)"
)

### membership is sparse and non-exclusive -------------------------------------

ica_modules_dt <- get_modules(results)

expect_true(
  current = data.table::uniqueN(ica_modules_dt$gene) < ncol(data),
  info = "ICA class - membership is sparse, not every feature is assigned"
)

expect_true(
  current = ica_modules_dt[, .N, by = gene][N > 1, .N] > 0,
  info = "ICA class - a gene can load on several components and belong to several modules"
)

expect_true(
  current = all(c("loading", "sign", "z") %in% names(ica_modules_dt)),
  info = "ICA class - membership reports loading, sign and the threshold statistic"
)

expect_true(
  current = nrow(get_modules(get_results(ica_stabilised_results(
    ica_test,
    no_comp = 3L,
    ica_type = "logcosh",
    membership_params = params_module_membership(method = "fdr")
  )))) >
    0L,
  info = "ICA class - the fdr membership method also returns modules"
)

# The fixtures keep the raw ICA `S` orientation (k x genes), hence the transpose
# on the way in and on the way out.
if (identical(Sys.getenv("REGEN"), "1")) {
  qs2::qs_save(
    t(get_factors(results, which = "gene_loadings")),
    "./test_data/ica_s_matrix_logcosh.qs"
  )
  qs2::qs_save(
    get_factors(results, which = "sample_activity"),
    "./test_data/ica_a_matrix_logcosh.qs"
  )
}

### class (exp) ----------------------------------------------------------------

ica_test <- ica_stabilised_results(ica_test, no_comp = 3L, ica_type = "exp")

results <- get_results(ica_test)

expect_equal(
  current = get_factors(results, which = "gene_loadings"),
  target = t(expected_ica_s_exp_mat),
  info = paste("ICA class - gene_loadings (exp, transposed S)")
)

expect_equal(
  current = get_factors(results, which = "sample_activity"),
  target = expected_ica_a_exp_mat,
  info = paste("ICA class - sample_activity (exp, A)")
)

if (identical(Sys.getenv("REGEN"), "1")) {
  qs2::qs_save(
    t(get_factors(results, which = "gene_loadings")),
    "./test_data/ica_s_matrix_exp.qs"
  )
  qs2::qs_save(
    get_factors(results, which = "sample_activity"),
    "./test_data/ica_a_matrix_exp.qs"
  )
}
