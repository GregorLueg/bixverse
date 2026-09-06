# Calculate the Hedge's G effect between two matrices

This function takes two matrices in and calculate on a per column basis
the Hedge's G effect size and the standard error. These results can be
subsequently used for meta-analyses or other approaches.

## Usage

``` r
calculate_effect_size(
  mat_a,
  mat_b,
  small_sample_correction = NULL,
  .verbose = TRUE
)
```

## Arguments

- mat_a:

  Numerical matrix. Contains the values for group a. Assumes that rows =
  samples, and columns = features.

- mat_b:

  Numerical matrix. Contains the values for group b.

- small_sample_correction:

  Can be NULL (automatic determination if a small sample size correction
  should be applied) or Boolean.

- .verbose:

  Boolean that controls verbosity of the function.

## Value

x, robustly scaled.

## Examples

``` r
# Hedge's G between two groups of ten samples
set.seed(42)
mat_a <- matrix(rnorm(100), nrow = 10, ncol = 10)
mat_b <- matrix(rnorm(100, mean = 1), nrow = 10, ncol = 10)
colnames(mat_a) <- colnames(mat_b) <- sprintf("gene_%i", 1:10)
res <- calculate_effect_size(mat_a, mat_b, .verbose = FALSE)
head(res$effect_sizes)
#> [1] -0.7178191 -0.6166706 -0.3550481 -1.4490642 -0.9636552 -0.6221444
```
