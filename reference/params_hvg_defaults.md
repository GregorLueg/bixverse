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

- loess_span - In case of

- search_budget - Search budget per tree for Annoy. Defaults to `100L`.

- n_trees - Number of trees to generate for Annoy. Defaults to `100L`.
