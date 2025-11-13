# Run PCA for single cell

This function will run PCA (option of full SVD and randomised SVD for
now) on the detected highly variable genes.

## Usage

``` r
calculate_pca_sc(
  object,
  no_pcs,
  randomised_svd = TRUE,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- no_pcs:

  Integer. Number of PCs to calculate.

- randomised_svd:

  Boolean. Shall randomised SVD be used. Faster, but less precise.

- seed:

  Integer. Controls reproducibility. Only relevant if
  `randomised_svd = TRUE`.

- .verbose:

  Boolean. Controls verbosity and returns run times.

## Value

The function will add the PCA factors, loadings and singular values to
the object cache in memory.
