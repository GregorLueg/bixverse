# gsva tests -------------------------------------------------------------------

## gaussian kernel -------------------------------------------------------------

### data -----------------------------------------------------------------------

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

expected_gaussian <- qs2::qs_read("./test_data/gsva_gaussian_res.qs")

### tests ----------------------------------------------------------------------

gaussian <- rs_gsva(
  exp = data,
  gs_list = random_gene_sets,
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

# this should return a matrix of NAs

gaussian_v2 <- rs_gsva(
  exp = data,
  gs_list = random_gene_sets_2,
  tau = 1.0,
  gaussian = TRUE,
  max_diff = TRUE,
  abs_rank = FALSE,
  timings = FALSE
)

expect_true(
  current = all(is.na(gaussian_v2)),
  info = "GSVA Gaussian kernel test - No overlapping genes"
)

### comparison against GSVA ----------------------------------------------------

if (requireNamespace("GSVA", quietly = TRUE)) {
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
}

## poisson kernel --------------------------------------------------------------

### data -----------------------------------------------------------------------

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

expected_poisson <- qs2::qs_read("./test_data/gsva_poisson_res.qs")

### tests ----------------------------------------------------------------------

poisson <- rs_gsva(
  exp = data_2,
  gs_list = random_gene_sets,
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
  gs_list = random_gene_sets_2,
  tau = 1.0,
  gaussian = TRUE,
  max_diff = TRUE,
  abs_rank = FALSE,
  timings = FALSE
)

expect_true(
  current = all(is.na(poisson_v2)),
  info = "GSVA Poisson kernel test - No overlapping genes"
)

### comparison against GSVA ----------------------------------------------------

if (requireNamespace("GSVA", quietly = TRUE)) {
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
}
