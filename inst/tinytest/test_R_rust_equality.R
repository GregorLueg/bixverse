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

cor_class <- bixverse:::upper_triangular_cor_mat$new(
  cor_coef = cor_data,
  features = colnames(mat),
  shift = 1L
)

expect_equal(
  current = cor_class$get_cor_matrix(.verbose = FALSE),
  target = cor(mat),
  info = "Upper triangle class test Rust <> R"
)
