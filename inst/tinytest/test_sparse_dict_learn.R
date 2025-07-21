# sparse dict learning tests ---------------------------------------------------

## data ------------------------------------------------------------------------

### simple test data -----------------------------------------------------------

set.seed(42L)

# parameters
n_samples <- 15
n_features <- 30
true_dict_size <- 4

true_dict <- matrix(0, nrow = n_samples, ncol = true_dict_size)

# pattern 1: constant
true_dict[, 1] <- 1.0

# pattern 2: linear
true_dict[, 2] <- seq(0, 1, length.out = n_samples)

# pattern 3: step function
true_dict[, 3] <- c(
  rep(1, n_samples %/% 2),
  rep(-1, n_samples - n_samples %/% 2)
)

# pattern 4: alternating
true_dict[, 4] <- rep(c(1, -1), length.out = n_samples)

# Normalise columns
for (j in 1:true_dict_size) {
  norm_val <- sqrt(sum(true_dict[, j]^2))
  if (norm_val > 1e-12) {
    true_dict[, j] <- true_dict[, j] / norm_val
  }
}

# Generate sparse coefficients
sparsity <- 2
true_coeffs <- matrix(0, nrow = true_dict_size, ncol = n_features)

for (j in 1:n_features) {
  # Randomly select which atoms to use
  active_atoms <- sample(1:true_dict_size, sparsity)
  true_coeffs[active_atoms, j] <- runif(sparsity, 0.5, 2.0)
}

# Generate synthetic data: Y = D * X + noise
synthetic_data <- true_dict %*% true_coeffs
noise <- matrix(
  runif(n_samples * n_features, -0.05, 0.05),
  nrow = n_samples,
  ncol = n_features
)
synthetic_data <- synthetic_data + noise

### synthetic data resembling biology ------------------------------------------

synthetic_data_2 <- generate_gene_module_data()

# expected return data
expected_dictionary <- qs2::qs_read("./test_data/dgrdl_dictionary.qs")
expected_coefficients <- qs2::qs_read("./test_data/dgrdl_coefficients.qs")
expected_feat_laplacian <- qs2::qs_read("./test_data/dgrdl_feat_laplacian.qs")
expected_sample_laplacian <- qs2::qs_read(
  "./test_data/dgrdl_sample_laplacian.qs"
)

## tests -----------------------------------------------------------------------

### simple data ----------------------------------------------------------------

res_simple <- rs_sparse_dict_dgrdl(
  x = synthetic_data,
  dgrdl_params = params_dgrdl(
    sparsity = 2L,
    dict_size = 4L,
    alpha = 0.01,
    beta = 0.01,
    max_iter = 8L,
    k_neighbours = 3L,
    admm_iter = 5L,
    rho = 1.0
  ),
  verbose = FALSE,
  seed = 10101L
)

reconstruction_simple <- res_simple$dictionary %*% res_simple$coefficients

reconstruction_error_simple <- norm(
  synthetic_data - reconstruction_simple,
  "F"
) /
  norm(synthetic_data, "F")

expect_equal(
  current = dim(reconstruction_simple),
  target = dim(synthetic_data),
  info = "DGRDL synthetic data 1 dimensions"
)

expect_true(
  current = reconstruction_error_simple < 0.2,
  info = "DGRDL synthetic data 1 reconstruction error"
)

avg_sparsity <- mean(colSums(abs(res_simple$coefficients) > 1e-6))

expect_equal(
  current = avg_sparsity,
  target = 2,
  info = "DGRDL synthetic data 1 sparsity"
)

### 'biological' data ----------------------------------------------------------

res_bio <- rs_sparse_dict_dgrdl(
  x = synthetic_data_2,
  dgrdl_params = params_dgrdl(
    sparsity = 3L,
    dict_size = 8L,
    alpha = 0.3,
    beta = 0.5,
    max_iter = 10L,
    k_neighbours = 3L,
    admm_iter = 5L,
    rho = 1.0
  ),
  verbose = FALSE,
  seed = 10101L
)

reconstruction_bio <- res_bio$dictionary %*% res_bio$coefficients
reconstruction_error_bio <- norm(synthetic_data_2 - reconstruction_bio, "F") /
  norm(synthetic_data_2, "F")

expect_true(
  current = reconstruction_error_bio < 0.2,
  info = "DGRDL synthetic data 2 reconstruction error"
)

expect_equal(
  current = res_bio$dictionary,
  target = expected_dictionary,
  info = "DGRDL synthetic data 2 expected dictionary"
)

expect_equal(
  current = res_bio$coefficients,
  target = expected_coefficients,
  info = "DGRDL synthetic data 2 expected coefficients"
)

expect_equal(
  current = res_bio$feature_laplacian,
  target = expected_feat_laplacian,
  info = "DGRDL synthetic data 2 expected feature laplacian"
)

expect_equal(
  current = res_bio$sample_laplacian,
  target = expected_sample_laplacian,
  info = "DGRDL synthetic data 2 expected sample laplacian"
)
