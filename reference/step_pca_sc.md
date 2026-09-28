# Pipeline step: PCA

Wraps
[`calculate_pca_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_pca_sc.md)
as an `ScStep`.

## Usage

``` r
step_pca_sc(
  no_pcs = 30L,
  pca_params = params_sc_pca(),
  sparse_svd = NULL,
  hvg = NULL,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- no_pcs:

  Integer. Number of PCs to calculate.

- pca_params:

  Named list. Controls the parameters to be used for the PCA calculation
  which is single cell-specific, see
  [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md)

- sparse_svd:

  Optional boolean. Solve on the sparse data with implicit centring and
  scaling instead of materialising the dense scaled matrix. `NULL` picks
  the sparse path, which is faster and lighter on memory for every
  solver, or the dense one for `residuals = TRUE`. `FALSE` uses the
  dense path, but only below 500,000 cells. Not used for `MetaCells`.

- hvg:

  Optional integer. If you want to provide your own HVG genes.
  Otherwise, the function will default to what is found in
  [`get_hvg()`](https://gregorlueg.github.io/bixverse/reference/get_hvg.md).
  Please provide 1-indexed genes here! If you provide these, the
  internal HVG will be overwritten.

- seed:

  Integer. Controls reproducibility. Only relevant for
  `svd_solver = "randomised"` and the Lanczos start vector.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

An `ScStep`.

## Examples

``` r
# PCA restricted to whatever the HVG step selected
step_hvg_sc(hvg_no = 30L) %>>% step_pca_sc(no_pcs = 10L)
#> <ScPipeline> 2 steps
#>   1. hvg  hvg_no = 30L, hvg_params = <list>, streaming = NULL, .verbose = TRUE
#>   2. pca  no_pcs = 10L, pca_params = <list>, sparse_svd = NULL, hvg = NULL, seed = 42L, .verbose = TRUE
```
