# PCA on MetaCells (sparse data)

**\[experimental\]** Calculates PCA for MetaCells or more generally
speaking sparse data. This is happening in-memory compared to the
(usually much) larger single cell data sets. The matrix is densified,
optionally CLR transformed and scaled according to `pca_params` before
the SVD.

## Usage

``` r
rs_mc_pca(sparse_data, no_pcs, pca_params, clr_offsets, seed, verbose)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `cs_type`. Shape is (metacells, genes), holding the
  normalised counts of the genes to use.

- no_pcs:

  Integer. Number of PCs to return.

- pca_params:

  Named list. Contains the parameters to use for this PCA run, see
  [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md).

- clr_offsets:

  Optional numeric. One offset per meta cell for the `PFlogPF`
  normalisation from Booeshaghi, et al., computed against the full gene
  panel. Required if `pca_params$clr` is `TRUE`, ignored otherwise.

- seed:

  Integer. Random seed for the randomised SVD.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- scores - The samples projected on the PCA space (solved via sparse
  SVD).

- loadings - The loadings of the features for the PCA (solved via sparse
  SVD).

- singular_values - The singular values for the PCA (solved via sparse
  SVD).

## References

Booeshaghi, et al., bioRxive, 2026.
