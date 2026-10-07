# Wrapper function for HVG detection parameters.

Wrapper function for HVG detection parameters.

## Usage

``` r
params_sc_hvg(
  method = c("vst", "meanvarbin", "dispersion", "scran", "residual"),
  loess_span = 0.3,
  num_bin = 20L,
  bin_method = c("equal_width", "equal_freq"),
  scran = list()
)
```

## Arguments

- method:

  String. `"scran"` fits a weighted lowess trend to the variance of the
  log-expression against its mean and ranks genes by the residual, as in
  scran's `modelGeneVar()`. `"residual"` ranks genes by the residual
  variance of a model fitted with
  [`fit_residuals_sc()`](https://gregorlueg.github.io/bixverse/reference/fit_residuals_sc.md),
  and needs that fit on the object first. It also treats `hvg_no` as a
  per-group count and returns the union across groups, so a grouped fit
  can select more than `hvg_no` genes. One of
  `c("vst", "meanvarbin", "dispersion", "scran", "residual")`. Defaults
  to `"vst"`.

- loess_span:

  Numeric. The span parameter for the loess function that is used to
  standardise the variance for `method = "vst"`, and the lowess span of
  the trend for `method = "scran"`. Defaults to `0.3`.

- num_bin:

  Integer. Not yet implemented. Defaults to `20L`.

- bin_method:

  String. The binning method. One of `c("equal_width", "equal_freq")`.
  Defaults to `"equal_width"`.

- scran:

  List. Optional overrides for the `method = "scran"` trend. See
  [`params_hvg_scran_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_hvg_scran_defaults.md)
  for available parameters: `mean_filter`, `min_mean`, `transform`,
  `use_min_width`, `min_width` and `min_window_count`. See
  [`params_hvg_scran_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_hvg_scran_defaults.md)
  for the available elements. Defaults to
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A named list with the following elements:

- method - String. `"scran"` fits a weighted lowess trend to the
  variance of the log-expression against its mean and ranks genes by the
  residual, as in scran's `modelGeneVar()`. `"residual"` ranks genes by
  the residual variance of a model fitted with
  [`fit_residuals_sc()`](https://gregorlueg.github.io/bixverse/reference/fit_residuals_sc.md),
  and needs that fit on the object first. It also treats `hvg_no` as a
  per-group count and returns the union across groups, so a grouped fit
  can select more than `hvg_no` genes. One of
  `c("vst", "meanvarbin", "dispersion", "scran", "residual")`. Defaults
  to `"vst"`.

- loess_span - Numeric. The span parameter for the loess function that
  is used to standardise the variance for `method = "vst"`, and the
  lowess span of the trend for `method = "scran"`. Defaults to `0.3`.

- num_bin - Integer. Not yet implemented. Defaults to `20L`.

- bin_method - String. The binning method. One of
  `c("equal_width", "equal_freq")`. Defaults to `"equal_width"`.

- The elements of
  [`params_hvg_scran_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_hvg_scran_defaults.md),
  overridden by `scran`, spliced in at this position.
