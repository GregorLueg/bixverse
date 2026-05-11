# Helper function to generate HVG defaults

Helper function to generate HVG defaults

## Usage

``` r
params_hvg_defaults()
```

## Value

A list with default parameters for kNN searches. Following parameters:

- min_gene_var_pctl - Which percentile of the highly variable genes to
  include. Defaults to `0.7`.

- hvg_method - Which method to use to identify HVG. Defaults to `"vst"`.

- loess_span - In case of `"vst"` the span of the loess function.

- clip_max - The maximum clipping value (optional).

- n_bins - The number of bins to use for the `"meanvarbin"` HVG
  detection.

- binning_strategy - Which binning strategy to use for `"meanvarbin"`.
