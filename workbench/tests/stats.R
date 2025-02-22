rextendr::document()
devtools::document()
devtools::check()

samples = 100
features = 10000

?diffuse_seed_nodes

.calculate_var <- function(es, grp1n, grp2n) {
  # Checks
  checkmate::qassert(es, "n")
  checkmate::qassert(grp1n, "i")
  checkmate::qassert(grp2n, "i")
  # Function body
  (grp1n + grp2n) / (grp1n * grp2n) + (es * es) / (2 * (grp1n + grp2n))
}



set.seed(123)
x_a <- matrix(rnorm(samples * features, 3, 1), nrow = samples, ncol = features)
x_b <- matrix(rnorm(samples * features, 1, 1), , nrow = samples, ncol = features)
colnames(x_a) <- colnames(x_b) <- sprintf("feature_%i", 1:features)


x_a[1:5, 1:5]
x_b[1:5, 1:5]

tictoc::tic()
grp1m = colMeans(x_a)
grp1sd = matrixStats::colSds(x_a)
grp1n = nrow(x_a)
grp2m = colMeans(x_b)
grp2sd = matrixStats::colSds(x_b)
grp2n = nrow(x_b)

totaln <- grp1n + grp2n

dm <- grp1m - grp2m

sd_pooled <- sqrt((grp1sd ^ 2 * (grp1n - 1) + grp2sd ^ 2 *
                     (grp2n - 1)) / (grp1n + grp2n - 2))

es <- dm / sd_pooled
var <- .calculate_var(es, grp1n, grp2n)
se <- sqrt(var)
tictoc::toc()

?rs_hedges_g

tictoc::tic()
rust_results <- calculate_effect_size(x_a, x_b)
tictoc::toc()

rust_results$effect_sizes

rust_results$standard_errors[1:5]

se[1:5]

plot(
  es,
  rust_results$effect_sizes
)

es[1:10]
rust_results$effect_sizes[1:10]

