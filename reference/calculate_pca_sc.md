# Run PCA for single cell

This function will run PCA (option of full SVD and randomised SVD for
now) on the detected highly variable genes.

## Usage

``` r
calculate_pca_sc(
  object,
  no_pcs,
  randomised_svd = TRUE,
  sparse_svd = FALSE,
  hvg = NULL,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- no_pcs:

  Integer. Number of PCs to calculate.

- randomised_svd:

  Boolean. Shall randomised SVD be used. Faster, but less precise.

- sparse_svd:

  Boolean. Shall sparse solvers be used that do not do scaling. If set
  to yes, in the case of `random_svd = FALSE`, Lanczos iterations are
  used to solve the sparse SVD. With `random_svd = TRUE`, the sparse
  initial matrix is multiplied with the random matrix, yielding a much
  smaller dense matrix that does not increase the memory pressure
  massively.

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

  Boolean. Controls verbosity and returns run times.

## Value

The function will add the PCA factors, loadings and singular values to
the object cache in memory.
