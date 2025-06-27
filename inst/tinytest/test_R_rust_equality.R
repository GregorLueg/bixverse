# correlation, covariance, equality tests --------------------------------------

## simple versions -------------------------------------------------------------

set.seed(123)
mat <- matrix(data = rnorm(100), nrow = 10, ncol = 10)
rownames(mat) <- sprintf("sample_%i", 1:10)
colnames(mat) <- sprintf("feature_%i", 1:10)

# Pearson
expect_equivalent(
  current = rs_cor(mat, spearman = FALSE),
  target = cor(mat),
  info = "Correlation equivalence test Rust <> R"
)
# Spearman
expect_equivalent(
  current = rs_cor(mat, spearman = TRUE),
  target = cor(mat, method = "spearman"),
  info = "Spearman Correlation equivalence test Rust <> R"
)
# Co-variance
expect_equivalent(
  current = rs_covariance(mat),
  target = cov(mat),
  info = "Covariance equivalence test Rust <> R"
)

## upper triangle versions -----------------------------------------------------

# Check if the upper triangle class behaves as expected
cor_data <- rs_cor_upper_triangle(mat, spearman = FALSE, shift = 1L)

cor_class <- bixverse:::upper_triangular_sym_mat$new(
  values = cor_data,
  features = colnames(mat),
  shift = 1L
)

expect_equal(
  current = cor_class$get_sym_matrix(.verbose = FALSE),
  target = cor(mat),
  info = "Upper triangle class test Rust <> R"
)

## covariance to cor -----------------------------------------------------------

cov_data <- rs_covariance(mat)

expect_equivalent(
  current = rs_cov2cor(cov_data),
  target = cor(mat),
  info = "cov2cor equivalence test Rust <> R"
)


# hypergeom distributions ------------------------------------------------------

m <- 10
n <- 7
k <- 8
x <- 0:(k + 1)

rust_vals <- purrr::map_dbl(
  x,
  ~ {
    rs_phyper(
      q = .x,
      m = m,
      n = n,
      k = k
    )
  }
)

r_vals <- purrr::map_dbl(
  x,
  ~ {
    phyper(
      q = .x,
      m = m,
      n = n,
      k = k,
      lower.tail = F
    )
  }
)

expect_equal(
  current = rust_vals,
  target = r_vals,
  info = "Hypergeometric test values for Rust <> R."
)

# pca --------------------------------------------------------------------------

r_pca_res <- prcomp(mat)

rs_pca_res <- rs_prcomp(mat, scale = FALSE)

# Absolute as signs can flip
expect_equivalent(
  current = abs(rs_pca_res$scores),
  target = abs(r_pca_res$x),
  info = "Factors for PCA (via SVD) for Rust <> R."
)

expect_equivalent(
  current = abs(rs_pca_res$v),
  target = abs(r_pca_res$rotation),
  info = "Loadings/rotation for PCA (via SVD) for Rust <> R."
)

expect_equal(
  current = rs_pca_res$s,
  target = r_pca_res$sdev,
  info = "Standard deviation of the eigenvalues for PCA (via SVD) for Rust <> R."
)
