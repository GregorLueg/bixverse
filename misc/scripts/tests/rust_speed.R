devtools::load_all()

rextendr::document()

# Improvements in correlation speeds

set.seed(123L)
no_samples <- 1000
no_features <- 1000
random_data <- matrix(
  rnorm(no_samples * no_features),
  nrow = no_samples,
  ncol = no_features
)

# Co-variance
microbenchmark::microbenchmark(
  cov_r <- cov(random_data),
  cov_rust <- rs_covariance(random_data),
  times = 20L
)

# Pearson correlation
microbenchmark::microbenchmark(
  cor_res_r <- cor(random_data),
  cor_res_rust <- rs_cor(random_data, spearman = FALSE),
  times = 20L
)


# Spearman correlation
microbenchmark::microbenchmark(
  cor_res_r <- cor(random_data, method = "spearman"),
  cor_res_rust <- rs_cor(random_data, spearman = TRUE),
  times = 20L
)


# Larger data matrix
set.seed(123L)
no_samples <- 1000
no_features <- 5000
random_data_2 <- matrix(
  rnorm(no_samples * no_features),
  nrow = no_samples,
  ncol = no_features
)

# Pearson correlation
microbenchmark::microbenchmark(
  cor_res_r_2 <- cor(random_data_2),
  cor_res_rust_2 <- rs_cor(random_data_2, spearman = FALSE),
  times = 20L
)
