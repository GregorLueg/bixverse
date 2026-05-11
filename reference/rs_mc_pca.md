# PCA on MetaCells (sparse data)

Calculates PCA for MetaCells or more generally speaking sparse data.
This is happening in-memory compared to the (usually much) larger single
cell data sets.

## Usage

``` r
rs_mc_pca(sparse_data, no_pcs, random_svd, seed)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `format`.

- no_pcs:

  Integer. Number of PCs to return.

- random_svd:

  Boolean. Shall randomised SVD be used.

- seed:

  Integer. Random seed for the randomised SVD.

## Value

A list with with the following items

- scores - The samples projected on the PCA space (solved via sparse
  SVD).

- loadings - The loadings of the features for the PCA (solved via sparse
  SVD).

- singular_values - The singular values for the PCA (solved via sparse
  SVD).
