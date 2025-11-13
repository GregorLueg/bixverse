# Calculate the Hedge G effect between two matrices

This function takes two matrices in and calculate on a per column basis
the Hedge's G effect size and the standard error. These results can be
subsequently used for meta-analyses or other approaches.

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

x, robustly scaled.
