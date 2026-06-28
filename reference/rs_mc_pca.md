# PCA on MetaCells (sparse data)

**\[experimental\]** Calculates PCA for MetaCells or more generally
speaking sparse data. This is happening in-memory compared to the
(usually much) larger single cell data sets.

## Usage

``` r
rs_mc_pca(sparse_data, no_pcs, pca_params, clr_offsets, seed)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `format`.

- no_pcs:

  Integer. Number of PCs to return.

- pca_params:

  Named list. Contains the parameters to use for this PCA run.

- clr_offsets:

  Optional numeric. If you wish to use the `PFlogPF` normalisation prior
  to PCA from Booeshaghi, et al.

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

## References

Booeshaghi, et al., bioRxive, 2026.
