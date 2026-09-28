# Run PCA for single cell

This function will run PCA on the detected highly variable genes. The
solver is set via `svd_solver` in
[`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md);
the default builds the gene x gene cross-product from the sparse data,
which is exact and the fastest option for a few thousand HVGs.

## Usage

``` r
calculate_pca_sc(
  object,
  no_pcs,
  pca_params = params_sc_pca(),
  sparse_svd = NULL,
  hvg = NULL,
  seed = 42L,
  residuals = FALSE,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` (or potentially other) class.

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

- residuals:

  Boolean. Run the PCA on the Pearson residuals of a model fitted with
  [`fit_residuals_sc()`](https://gregorlueg.github.io/bixverse/reference/fit_residuals_sc.md)
  instead of the stored normalised layer. Needs
  `params_sc_pca(normalise_variance = FALSE, clr = FALSE)`, since the
  residuals already carry the signal as variance, and does not support
  `sparse_svd`: a residual column is dense even where the counts are
  not. On this dense path the default `svd_solver = "covariance"` pays
  for a full gene x gene cross-product and is slower than
  `svd_solver = "randomised"`, so pick the latter in
  [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md)
  when speed matters. Not supported for `MetaCells`.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

The function will add the PCA factors, loadings and singular values to
the object cache in memory.

## Examples

``` r
# PCA on the highly variable genes
sc <- demo_single_cells(prepped = FALSE)
sc <- find_hvg_sc(sc, hvg_no = 30L, .verbose = FALSE)
sc <- calculate_pca_sc(sc, no_pcs = 10L, .verbose = FALSE)
dim(get_pca_factors(sc))
#> [1] 500  10

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
