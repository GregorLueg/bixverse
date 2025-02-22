rextendr::document()
devtools::document()
devtools::check()

samples <- 100
features <- 10000

set.seed(123)
x_a <- matrix(rnorm(samples * features, 3, 1), nrow = samples, ncol = features)
x_b <- matrix(rnorm(samples * features, 1, 1), , nrow = samples, ncol = features)
colnames(x_a) <- colnames(x_b) <- sprintf("feature_%i", 1:features)

tictoc::tic()
rust_results <- calculate_effect_size(x_a, x_b)
tictoc::toc()
