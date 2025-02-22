library(devtools)

devtools::document()
devtools::load_all()
rextendr::document()


syn_data = synthetic_signal_matrix()

X = t(syn_data$mat)

X[1:5, 1:5]

cor(X)[1:5, 1:5]

rs_cor(X, spearman = FALSE)[1:5, 1:5]

X_a <- matrix(data = rnorm(1000 * 1000, 1, 2), nrow = 1000, ncol = 1000)
X_b <- matrix(data = rnorm(1000 * 1000), nrow = 1000, ncol = 1000)

diff_cor_res <- rs_differential_cor(X_a, X_b, spearman = TRUE)

tictoc::tic()
cor_X <- cor(X_a, method = "spearman")
tictoc::toc()

cor_X[1:5, 1:5]

tictoc::tic()
cor_X_rs <- rs_cor(X_a, spearman = TRUE)
tictoc::toc()

cor_X_rs[1:5, 1:5]

n <- nrow(X_scaled)

cor_m <- t(X_scaled) %*% X_scaled / (n - 1)

cor <- (m1 - row_means) / row_stds

cor_m[1:5, 1:5]
