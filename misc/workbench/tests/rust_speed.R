devtools::load_all()


# Improvements in correlation speeds

set.seed(123L)
no_samples <- 1000
no_features <- 1000
random_data <- matrix(rnorm(no_samples * no_features),
                      nrow = no_samples,
                      ncol = no_features)

# Pearson correlation
microbenchmark::microbenchmark(
  cor_res_r <- cor(random_data),
  cor_res_rust <- rs_cor(random_data, spearman = FALSE),
  times = 10L
)

# Co-variance
microbenchmark::microbenchmark(cov_r <- cov(random_data),
                               cov_rust <- rs_covariance(random_data),
                               times = 10L)

# Spearman correlation
microbenchmark::microbenchmark(
  cor_res_r <- cor(random_data, method = 'spearman'),
  cor_res_rust <- rs_cor(random_data, spearman = TRUE),
  times = 10L
)
