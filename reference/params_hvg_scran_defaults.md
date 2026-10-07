# Helper function to generate the scran HVG trend defaults

Trend parameters for `method = "scran"` in
[`params_sc_hvg()`](https://gregorlueg.github.io/bixverse/reference/params_sc_hvg.md).
They mirror scrapper's `fitVarianceTrend()` defaults. The lowess span
comes from `loess_span` in
[`params_sc_hvg()`](https://gregorlueg.github.io/bixverse/reference/params_sc_hvg.md).

## Usage

``` r
params_hvg_scran_defaults()
```

## Value

A named list with the following elements:

- mean_filter - Boolean. Shall genes below `min_mean` be left out of the
  trend fit. Defaults to `TRUE`.

- min_mean - Numeric. Minimum mean log-expression for a gene to enter
  the trend fit. Defaults to `0.1`.

- transform - Boolean. Shall the variances be fourth-root transformed
  before the fit. Defaults to `TRUE`.

- use_min_width - Boolean. Shall the lowess window be defined by
  `min_width` and `min_window_count` instead of the span. Defaults to
  `FALSE`.

- min_width - Numeric. Minimum window width, only used with
  `use_min_width = TRUE`. Defaults to `1.0`.

- min_window_count - Integer. Minimum number of genes per window, only
  used with `use_min_width = TRUE`. Defaults to `200L`.
