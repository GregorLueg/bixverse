# Reject PCA settings the residual path cannot honour

Both are refused rather than overridden. They are settings the caller
passed explicitly, and quietly changing them would make the recorded
parameters disagree with what was actually run.

## Usage

``` r
.assert_residual_pca_params(pca_params)
```

## Arguments

- pca_params:

  List. See
  [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md).

## Value

Invisibly `TRUE`; called for the error.
