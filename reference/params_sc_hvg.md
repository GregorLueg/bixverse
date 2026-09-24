# Wrapper function for HVG detection parameters.

Wrapper function for HVG detection parameters.

## Usage

``` r
params_sc_hvg(
  method = c("vst", "meanvarbin", "dispersion", "residual"),
  loess_span = 0.3,
  num_bin = 20L,
  bin_method = c("equal_width", "equal_freq")
)
```

## Arguments

- method:

  String. `"residual"` ranks genes by the residual variance of a model
  fitted with
  [`fit_residuals_sc()`](https://gregorlueg.github.io/bixverse/reference/fit_residuals_sc.md),
  and needs that fit on the object first. It also treats `hvg_no` as a
  per-group count and returns the union across groups, so a grouped fit
  can select more than `hvg_no` genes. One of
  `c("vst", "meanvarbin", "dispersion", "residual")`. Defaults to
  `"vst"`.

- loess_span:

  Numeric. The span parameter for the loess function that is used to
  standardise the variance for `method = "vst"`. Defaults to `0.3`.

- num_bin:

  Integer. Not yet implemented. Defaults to `20L`.

- bin_method:

  String. The binning method. One of `c("equal_width", "equal_freq")`.
  Defaults to `"equal_width"`.

## Value

A named list with the following elements:

- method - String. `"residual"` ranks genes by the residual variance of
  a model fitted with
  [`fit_residuals_sc()`](https://gregorlueg.github.io/bixverse/reference/fit_residuals_sc.md),
  and needs that fit on the object first. It also treats `hvg_no` as a
  per-group count and returns the union across groups, so a grouped fit
  can select more than `hvg_no` genes. One of
  `c("vst", "meanvarbin", "dispersion", "residual")`. Defaults to
  `"vst"`.

- loess_span - Numeric. The span parameter for the loess function that
  is used to standardise the variance for `method = "vst"`. Defaults to
  `0.3`.

- num_bin - Integer. Not yet implemented. Defaults to `20L`.

- bin_method - String. The binning method. One of
  `c("equal_width", "equal_freq")`. Defaults to `"equal_width"`.
