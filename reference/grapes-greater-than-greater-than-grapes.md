# Append a step to a pipeline

`%>>%` chains pipeline steps. Either side can be a `ScStep` or a
`ScPipeline`; the result is always a `ScPipeline`.

## Usage

``` r
lhs %>>% rhs
```

## Arguments

- lhs:

  `ScPipeline` or `ScStep`.

- rhs:

  `ScStep`.

## Value

A `ScPipeline`.

## Examples

``` r
# either side may be a step, the result is always a pipeline
step_hvg_sc(hvg_no = 30L) %>>% step_pca_sc(no_pcs = 10L)
#> <ScPipeline> 2 steps
#>   1. hvg  hvg_no = 30L, hvg_params = <list>, streaming = NULL, .verbose = TRUE
#>   2. pca  no_pcs = 10L, pca_params = <list>, sparse_svd = FALSE, hvg = NULL, seed = 42L, .verbose = TRUE
```
