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

- hvg_method - String. Which method to use to identify HVG. One of
  `c("vst", "mvb", "dispersion")`. Defaults to `"vst"`.

- loess_span - Numeric. In case of `"vst"` the span of the loess
  function. Defaults to `0.3`.

- clip_max - Numeric or `NULL`. The maximum clipping value (optional).
  Defaults to `NULL`.

- n_bins - Integer. The number of bins to use for the `"mvb"` HVG
  detection. Defaults to `20L`.

- binning_strategy - String. Which binning strategy to use for `"mvb"`.
  One of `c("equal_width", "equal_frequency")`. Defaults to
  `"equal_width"`.
