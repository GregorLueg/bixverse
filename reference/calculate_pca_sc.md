# Run PCA for single cell

This function will run PCA on the detected highly variable genes. You
can use randomised SVD for speed and there is an option for sparse SVD
for very large data sets to avoid memory pressure.

## Usage

``` r
calculate_pca_sc(
  object,
  no_pcs,
  pca_params = params_sc_pca(),
  sparse_svd = FALSE,
  hvg = NULL,
  seed = 42L,
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

  Boolean. Shall sparse solvers be used that do not do scaling. If set
  to yes, in the case of `random_svd = FALSE`, Lanczos iterations are
  used to solve the sparse SVD. With `random_svd = TRUE`, the sparse
  initial matrix is multiplied with the random matrix, yielding a much
  smaller dense matrix that does not increase the memory pressure
  massively. Not used for `MetaCells`.

- hvg:

  Optional integer. If you want to provide your own HVG genes.
  Otherwise, the function will default to what is found in
  [`get_hvg()`](https://gregorlueg.github.io/bixverse/reference/get_hvg.md).
  Please provide 1-indexed genes here! If you provide these, the
  internal HVG will be overwritten.

- seed:

  Integer. Controls reproducibility. Only relevant if
  `randomised_svd = TRUE`.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

The function will add the PCA factors, loadings and singular values to
the object cache in memory.
