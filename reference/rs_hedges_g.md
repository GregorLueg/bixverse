# Calculate the Hedge's G effect

Calculates the Hedge's G effect for two sets of matrices. The function
assumes that rows = samples and columns = features. WARNING! Incorrect
use can cause kernel crashes. Wrapper around the Rust functions with
type checks are provided in the package.

## Usage

``` r
rs_hedges_g(mat_a, mat_b, small_sample_correction)
```

## Arguments

- mat_a:

  The matrix of samples and features in grp A for which to calculate the
  Hedge's G effect.

- mat_b:

  The matrix of samples and features in grp B for which to calculate the
  Hedge's G effect.

- small_sample_correction:

  Shall the small sample correction be applied.

## Value

Returns the harmonic sum according to the OT calculation.
