# Run Seurat CCA integration

This function implements the canonical correlation analysis (CCA) anchor
integration from Stuart, et al. For each pair of batches a shared CCA
embedding is computed, mutual nearest neighbours in that embedding
become anchors, these are filtered in gene space, scored by shared
neighbours and finally used to apply a kernel-weighted correction on the
union PCA embedding. Batches are merged in the order of their pairwise
anchor counts.

This port deviates from Seurat in two places. It skips the per-gene
`ScaleData` step and works from per-cell standardised log-normalised HVG
expression, and it never materialises the `N1 x N2` cross-product (the
canonical correlations come from a matrix-free randomised SVD). The
correction runs on the embedding, not on full log-expression. In
practice the anchor structure comes out close to identical at a fraction
of the memory.

## Usage

``` r
seurat_cca_sc(
  object,
  batch_column,
  batch_hvg_genes,
  cca_params = params_sc_seurat_cca(),
  use_precomputed_pca = FALSE,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- batch_column:

  String. The column with the batch information in the obs data of the
  class.

- batch_hvg_genes:

  Integer vector. These are the highly variable genes, identified by a
  batch-aware method. Please refer to
  [`find_hvg_batch_aware_sc()`](https://gregorlueg.github.io/bixverse/reference/find_hvg_batch_aware_sc.md)
  for more details. These genes have to be 0-indexed!

- cca_params:

  A list, please see
  [`params_sc_seurat_cca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_seurat_cca.md).
  The list has the following parameters:

  - num_cc - Integer. Number of canonical correlation dimensions. The
    effective rank used is `max(num_cc, dims)`.

  - dims - Integer. Dimensions used for the anchor kNN queries and size
    of the returned embedding.

  - k_anchor - Integer. Neighbourhood size for the anchor search.

  - k_filter - Integer. Neighbourhood size for the gene-space filter.

  - k_score - Integer. Neighbourhood size for the anchor scoring.

  - k_weight - Integer. Neighbourhood size for the kernel weights.

  - n_top_features - Integer. Top-loading genes for the gene-space
    filter.

  - l2_norm - Logical. L2-normalise the CCA embedding per cell.

  - sd - Numeric. Bandwidth divisor of the Gaussian kernel.

  - knn - List of kNN parameters. See
    [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
    for available parameters and their defaults.

  - pca - List of PCA parameters, see
    [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md)
    for available parameters and their defaults.

- use_precomputed_pca:

  Boolean. Should the PCA in the object be used if found. If you decide
  to do this, make sure that you have run the PCA on the batch-aware HVG
  ideally. Note that CCA still needs the PCA loadings for the gene-space
  filter, so this saves less work than it does for fastMNN.

- seed:

  Integer. Random seed.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

The object with the added `"cca"` embedding.

## References

Stuart, et al., Cell, 2019
