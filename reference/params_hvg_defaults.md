# Helper function to generate HVG defaults

Helper function to generate HVG defaults

## Usage

``` r
params_hvg_defaults()
```

## Value

A named list with the following elements:

- min_gene_var_pctl - Numeric. Which percentile of the highly variable
  genes to include. Defaults to `0.7`.

- hvg_method - String. Which method to use to identify HVG. `"scran"`
  runs with the default trend parameters, see
  [`params_hvg_scran_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_hvg_scran_defaults.md).
  One of `c("vst", "meanvarbin", "dispersion", "scran")`. Defaults to
  `"vst"`.

- loess_span - Numeric. In case of `"vst"` the span of the loess
  function. Defaults to `0.3`.

- clip_max - Numeric or `NULL`. The maximum clipping value (optional).
  Defaults to `NULL`.

- n_bins - Integer. The number of bins to use for the `"meanvarbin"` and
  `"dispersion"` HVG detection. Defaults to `20L`.

- binning_strategy - String. Which binning strategy to use for
  `"meanvarbin"` and `"dispersion"`. One of
  `c("equal_width", "equal_frequency")`. Defaults to `"equal_width"`.
