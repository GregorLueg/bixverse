# Construct the fitted residual model side-car

Thin wrapper over the list
[`rs_sc_fit_residuals()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_fit_residuals.md)
returns, with the provenance R knows about and Rust does not.

`covariate_names` is kept separately even though the covariates carry
their own names: the fitted coefficients are matched to covariates by
position, so a reordered selection would be applied to the wrong column.
Rust checks the order on every use, and holding the expected order here
lets R say what went wrong first.

## Usage

``` r
new_sc_residual_fit(res, method, params, group_column = NULL)
```

## Arguments

- res:

  List. The result of
  [`rs_sc_fit_residuals()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_fit_residuals.md).

- method:

  String. The method that was fitted.

- params:

  List. The parameters used.

- group_column:

  String or `NULL`. The grouping column, if any.

## Value

The `ScResidualFit` object.
