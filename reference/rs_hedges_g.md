# Calculate the Hedge's G effect

**\[experimental\]** Calculates the Hedge's G effect for two sets of
matrices. The function assumes that rows = samples and columns =
features.

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

A list with the following items:

- effect_sizes - Hedge's G effect size per feature.

- standard_errors - Standard error of the effect size per feature.
