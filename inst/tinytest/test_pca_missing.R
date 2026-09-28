# pca with missing values ------------------------------------------------------

## data ------------------------------------------------------------------------

# rank-2 signal plus small noise, 6 of 60 entries missing (10%)
set.seed(123L)
z <- matrix(rnorm(12L * 2L), 12L, 2L)
w <- matrix(rnorm(5L * 2L), 5L, 2L)
x <- z %*% t(w) + matrix(rnorm(12L * 5L, sd = 0.1), 12L, 5L)
na_idx <- c(3L, 17L, 26L, 38L, 45L, 58L)
x[na_idx] <- NA

# align each column of `a` to the sign of the matching column of `b`
align_sign <- function(a, b) {
  a %*% diag(sign(colSums(a * b)), ncol(a))
}

## expected values -------------------------------------------------------------

# pcaMethods 2.0.0 on R 4.5, run on `x` above:
# pca(x, method = "ppca", nPcs = 2L, center = TRUE, seed = 42L)
# pca(x, method = "bpca", nPcs = 2L, center = TRUE)
# imputed = completeObs(.)[na_idx]; values via dput(signif(., 8))

ppca_imp <- c(
  -1.5090233,
  0.22360094,
  -0.1758079,
  0.081361943,
  -1.1799677,
  0.43728413
)
ppca_scores <- matrix(
  c(
    -2.5485641,
    -1.426758,
    3.3661807,
    -3.4314961,
    -1.4204019,
    6.0079436,
    -0.83242698,
    -2.7823856,
    -0.44608288,
    -1.5919047,
    3.5616014,
    1.309711,
    -0.13887972,
    -0.093604016,
    0.79792647,
    1.5996668,
    0.29997215,
    -0.48624461,
    0.885158,
    -1.7394304,
    -1.7805126,
    -0.63811876,
    0.12336466,
    -0.50125977
  ),
  nrow = 12L
)
ppca_loadings <- matrix(
  c(
    -0.49423521,
    -0.5795168,
    0.30676344,
    -0.17161435,
    -0.54436802,
    0.52231058,
    -0.52521886,
    0.17635003,
    0.6479496,
    -0.019969644
  ),
  nrow = 5L
)
ppca_r2 <- c(0.8929193, 0.99786752)

bpca_imp <- c(
  -1.5137413,
  0.22553874,
  -0.17439236,
  0.082699753,
  -1.1462973,
  0.43412761
)
bpca_scores <- matrix(
  c(
    -0.88926086,
    -0.49788549,
    1.1783705,
    -1.1877475,
    -0.49391138,
    2.0919739,
    -0.28545881,
    -0.97936962,
    -0.16646272,
    -0.5580005,
    1.2423706,
    0.45391436,
    -0.13176596,
    -0.089432006,
    0.8143226,
    1.6888754,
    0.31973041,
    -0.53918774,
    0.92887875,
    -1.8025916,
    -1.8408371,
    -0.65818984,
    0.11072259,
    -0.53071376
  ),
  nrow = 12L
)
bpca_loadings <- matrix(
  c(
    -1.4072972,
    -1.6663845,
    0.88040232,
    -0.48230625,
    -1.5581077,
    0.49940888,
    -0.48638711,
    0.16134262,
    0.6104526,
    -0.010245646
  ),
  nrow = 5L
)
bpca_r2 <- c(0.90330535, 0.99773693)

## ppca ------------------------------------------------------------------------

# the start comes from the Rust RNG, not R's, so agreement is only at
# convergence (tol = 1e-5)
ppca_res <- run_ppca(x, ppca_params = params_ppca(n_pcs = 2L), .verbose = FALSE)

expect_true(ppca_res$converged, info = "PPCA converged")
expect_equal(
  ppca_res$completed[na_idx],
  ppca_imp,
  tolerance = 1e-4,
  info = "PPCA imputed values vs pcaMethods"
)
expect_identical(
  ppca_res$completed[-na_idx],
  x[-na_idx],
  info = "PPCA leaves observed entries untouched"
)
expect_equal(
  unname(align_sign(ppca_res$scores, ppca_scores)),
  ppca_scores,
  tolerance = 1e-4,
  info = "PPCA scores vs pcaMethods"
)
expect_equal(
  unname(align_sign(ppca_res$loadings, ppca_loadings)),
  ppca_loadings,
  tolerance = 1e-4,
  info = "PPCA loadings vs pcaMethods"
)
expect_equal(
  ppca_res$r2_cum,
  ppca_r2,
  tolerance = 1e-5,
  info = "PPCA cumulative R2 vs pcaMethods"
)
expect_equal(
  colnames(ppca_res$scores),
  c("PC1", "PC2"),
  info = "PPCA components are named"
)

## bpca ------------------------------------------------------------------------

# deterministic SVD start, so this should match pcaMethods to rounding
bpca_res <- run_bpca(x, bpca_params = params_bpca(n_pcs = 2L), .verbose = FALSE)

expect_true(bpca_res$converged, info = "BPCA converged")
expect_equal(
  bpca_res$completed[na_idx],
  bpca_imp,
  tolerance = 1e-6,
  info = "BPCA imputed values vs pcaMethods"
)
expect_identical(
  bpca_res$completed[-na_idx],
  x[-na_idx],
  info = "BPCA leaves observed entries untouched"
)
expect_equal(
  unname(align_sign(bpca_res$scores, bpca_scores)),
  bpca_scores,
  tolerance = 1e-6,
  info = "BPCA scores vs pcaMethods"
)
expect_equal(
  unname(align_sign(bpca_res$loadings, bpca_loadings)),
  bpca_loadings,
  tolerance = 1e-6,
  info = "BPCA loadings vs pcaMethods"
)
expect_equal(
  bpca_res$r2_cum,
  bpca_r2,
  tolerance = 1e-6,
  info = "BPCA cumulative R2 vs pcaMethods"
)

## missingness encoding --------------------------------------------------------

# NaN and NA both mark a missing entry
x_nan <- x
x_nan[na_idx] <- NaN
bpca_res_nan <- run_bpca(
  x_nan,
  bpca_params = params_bpca(n_pcs = 2L),
  .verbose = FALSE
)
expect_equal(
  bpca_res_nan$completed,
  bpca_res$completed,
  info = "NaN is treated like NA"
)

## input validation ------------------------------------------------------------

expect_error(
  run_ppca(as.vector(x), .verbose = FALSE),
  pattern = "Must be of type 'matrix'"
)

x_inf <- x
x_inf[1L, 1L] <- Inf
expect_error(
  run_bpca(x_inf, .verbose = FALSE),
  pattern = "Must be finite"
)

x_empty_row <- x
x_empty_row[2L, ] <- NA
expect_error(
  run_ppca(x_empty_row, .verbose = FALSE),
  pattern = "Rows entirely missing: 2\\."
)

x_empty_col <- x
x_empty_col[, 3L] <- NA
expect_error(
  run_bpca(x_empty_col, .verbose = FALSE),
  pattern = "Columns entirely missing: 3\\."
)

expect_error(
  run_bpca(x, bpca_params = params_bpca(n_pcs = 6L), .verbose = FALSE),
  pattern = "at most 5 are possible"
)

expect_error(
  run_ppca(x, ppca_params = list(n_pcs = 2L), .verbose = FALSE),
  pattern = "Names must include"
)

expect_error(
  params_ppca(n_pcs = 0L),
  pattern = "All elements must be >= 1"
)

expect_error(
  params_bpca(tol = -1),
  pattern = "All elements must be > 0"
)
