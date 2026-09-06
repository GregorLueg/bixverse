# Construct an empty single cell pipeline

Linear container of `ScStep`s. Append steps with `%>>%` and execute with
[`apply_pipeline()`](https://gregorlueg.github.io/bixverse/reference/apply_pipeline.md).
Pipelines are inert until applied; steps can be inspected via
`pipeline$steps`.

## Usage

``` r
sc_pipeline()
```

## Value

An empty `ScPipeline` object.

## Examples

``` r
# an empty container, filled with `%>>%`
sc_pipeline() %>>%
  step_hvg_sc(hvg_no = 30L) %>>%
  step_pca_sc(no_pcs = 10L)
#> <ScPipeline> 2 steps
#>   1. hvg  hvg_no = 30L, hvg_params = <list>, streaming = NULL, .verbose = TRUE
#>   2. pca  no_pcs = 10L, pca_params = <list>, sparse_svd = FALSE, hvg = NULL, seed = 42L, .verbose = TRUE
```
