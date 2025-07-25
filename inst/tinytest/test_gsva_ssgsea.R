# gsva / ssgsea tests ----------------------------------------------------------

## data ------------------------------------------------------------------------

### gaussian -------------------------------------------------------------------

set.seed(42L)

no_genes <- 100L
no_samples <- 10L
no_gs <- 5L

data <- matrix(
  data = rnorm(no_genes * no_samples, sd = 2),
  nrow = no_genes,
  ncol = no_samples
)

rownames(data) <- sprintf("gene_%i", 1:no_genes)
colnames(data) <- sprintf("sample_%i", 1:no_samples)

random_gene_sets <- purrr::map(
  1:no_gs,
  ~ {
    sample(rownames(data), 5)
  }
)
names(random_gene_sets) <- sprintf("gs_%i", 1:no_gs)

random_gene_sets_2 <- purrr::map(
  1:no_gs,
  ~ {
    sample(LETTERS, 5)
  }
)
names(random_gene_sets_2) <- sprintf("gs_%i", 1:no_gs)

### poisson --------------------------------------------------------------------

set.seed(123L)

no_genes <- 100L
no_samples <- 10L
no_gs <- 5L

data_2 <- matrix(
  data = rpois(no_genes * no_samples, lambda = 10),
  nrow = no_genes,
  ncol = no_samples
)

rownames(data_2) <- sprintf("gene_%i", 1:no_genes)
colnames(data_2) <- sprintf("sample_%i", 1:no_samples)

storage.mode(data_2) <- "numeric"

### expected data --------------------------------------------------------------

expected_gaussian <- qs2::qs_read("./test_data/gsva_gaussian_res.qs")
expected_poisson <- qs2::qs_read("./test_data/gsva_poisson_res.qs")

expected_gs_indices <- list(
  gs_1 = c(3, 11, 26, 44, 65),
  gs_2 = c(27, 70, 76, 94, 99),
  gs_3 = c(9, 40, 47, 64, 89),
  gs_4 = c(46, 56, 64, 66, 70),
  gs_5 = c(1, 8, 13, 58, 84)
)

expected_bad_matrix <- matrix(data = rep(NA, 10), nrow = 1, ncol = 10)
storage.mode(expected_bad_matrix) <- "numeric"

## tests gsva ------------------------------------------------------------------

### gene set indexing ----------------------------------------------------------

gs_indices <- rs_prepare_gsva_gs(
  feature_names = rownames(data),
  pathway_list = random_gene_sets,
  min_size = 5L,
  max_size = 500L
)

gs_indices_2 <- rs_prepare_gsva_gs(
  feature_names = rownames(data),
  pathway_list = random_gene_sets_2,
  min_size = 5L,
  max_size = 500L
)

expect_equal(
  current = gs_indices,
  target = expected_gs_indices,
  info = "GSVA - actual GS indices"
)

expect_true(
  current = purrr::is_empty(gs_indices_2),
  info = "GSVA - GS indices if nothing matches"
)

### gaussian -------------------------------------------------------------------

#### rust ----------------------------------------------------------------------

# direct rust version
gaussian <- rs_gsva(
  exp = data,
  gs_list = gs_indices,
  tau = 1.0,
  gaussian = TRUE,
  max_diff = TRUE,
  abs_rank = FALSE,
  timings = FALSE
)

expect_equal(
  current = gaussian,
  target = expected_gaussian,
  info = "GSVA Gaussian kernel test"
)

# this should 1 row of NAs

gaussian_v2 <- rs_gsva(
  exp = data,
  gs_list = gs_indices_2,
  tau = 1.0,
  gaussian = TRUE,
  max_diff = TRUE,
  abs_rank = FALSE,
  timings = FALSE
)

expect_equal(
  current = gaussian_v2,
  target = expected_bad_matrix,
  info = "GSVA Gaussian kernel test - empty GS list"
)

#### r wrapper -----------------------------------------------------------------

r_gaussian <- calc_gsva(
  exp = data,
  pathways = random_gene_sets
)

expect_equivalent(
  current = r_gaussian,
  target = expected_gaussian,
  info = "GSVA Gaussian kernel test - R wrapper"
)

# check that errors are happening
expect_error(
  current = calc_gsva(
    exp = data,
    pathways = random_gene_sets,
    gsva_params = params_gsva(min_size = 1L)
  )
)

expect_error(
  current = calc_gsva(
    exp = data,
    pathways = random_gene_sets,
    gsva_params = params_gsva(max_size = 3L)
  )
)

### poisson --------------------------------------------------------------------

#### rust ----------------------------------------------------------------------

poisson <- rs_gsva(
  exp = data_2,
  gs_list = gs_indices,
  tau = 1.0,
  gaussian = FALSE,
  max_diff = TRUE,
  abs_rank = FALSE,
  timings = FALSE
)

expect_equal(
  current = poisson,
  target = expected_poisson,
  info = "GSVA Poisson kernel test"
)

# this should return a matrix of NAs

poisson_v2 <- rs_gsva(
  exp = data_2,
  gs_list = gs_indices_2,
  tau = 1.0,
  gaussian = TRUE,
  max_diff = TRUE,
  abs_rank = FALSE,
  timings = FALSE
)

expect_equal(
  current = poisson_v2,
  target = expected_bad_matrix,
  info = "GSVA Poisson kernel test - empty GS list"
)

#### r wrapper -----------------------------------------------------------------

r_poisson <- calc_gsva(
  exp = data_2,
  pathways = random_gene_sets,
  gaussian = FALSE
)

expect_equivalent(
  current = r_poisson,
  target = expected_poisson,
  info = "GSVA Poisson kernel test - R wrapper"
)


## tests ssgsea ----------------------------------------------------------------

### expected data --------------------------------------------------------------

expected_ssgsea <- qs2::qs_read("./test_data/ssgsea_res.qs")
expected_ssgesa_unnorm <- qs2::qs_read("./test_data/ssgsea_res_unnorm.qs")

### rust -----------------------------------------------------------------------

ssgsea_res <- rs_ssgsea(
  exp = data,
  gs_list = gs_indices,
  alpha = 0.25,
  normalise = TRUE,
  timings = FALSE
)

ssgsea_res_unnorm <- rs_ssgsea(
  exp = data,
  gs_list = gs_indices,
  alpha = 0.25,
  normalise = FALSE,
  timings = FALSE
)

expect_equal(
  current = ssgsea_res,
  target = expected_ssgsea,
  info = "ssGSEA Rust interface (normalised)"
)

expect_equal(
  current = ssgsea_res_unnorm,
  target = expected_ssgesa_unnorm,
  info = "ssGSEA Rust interface (unnormalised)"
)

### rust interface -------------------------------------------------------------

ssgsea_r_res <- calc_ssgsea(
  exp = data,
  pathways = random_gene_sets
)

expect_equivalent(
  current = ssgsea_r_res,
  target = expected_ssgsea,
  info = "ssGSEA - R wrapper"
)

## comparison against GSVA -----------------------------------------------------

if (requireNamespace("GSVA", quietly = TRUE)) {
  ### Gaussian version

  gsvaPar <- GSVA::gsvaParam(data, random_gene_sets)
  gsva_results <- as.matrix(GSVA::gsva(
    gsvaPar,
    verbose = FALSE
  ))

  # There are tiny differences due to Gaussian optimisation
  expect_equivalent(
    current = gaussian,
    target = gsva_results,
    info = "GSVA Rust vs R original implementation (Gaussian)",
    tolerance = 1e-3
  )

  # However, the correlations are super high, so okay
  expect_true(
    current = all(diag(cor(gaussian, gsva_results)) > 0.99),
    info = "GSVA Rust vs R original implementation correlation (Gaussian)"
  )

  ### Poisson version

  gsvaPar <- GSVA::gsvaParam(data_2, random_gene_sets, kcdf = "Poisson")
  gsva_results <- as.matrix(GSVA::gsva(
    gsvaPar,
    verbose = FALSE
  ))

  expect_equivalent(
    current = poisson,
    target = gsva_results,
    info = "GSVA Rust vs R original implementation (Poisson)"
  )

  # However, the correlations are super high, so okay
  expect_true(
    current = all(diag(cor(poisson, gsva_results)) > 0.99),
    info = "GSVA Rust vs R original implementation correlation (Poisson)"
  )

  ### ssGSEA

  ssgseaPar <- GSVA::ssgseaParam(data, random_gene_sets)
  ssgsea_res <- as.matrix(GSVA::gsva(
    ssgseaPar,
    verbose = FALSE
  ))

  expect_equivalent(
    current = ssgsea_r_res,
    target = ssgsea_res,
    info = "ssGSEA Rust vs R original implementation"
  )

  expect_true(
    current = all(diag(cor(ssgsea_r_res, ssgsea_res)) > 0.99),
    info = "ssGSEA Rust vs R original implementation correlation"
  )
}
