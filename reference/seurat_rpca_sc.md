# Run Seurat rPCA integration

This function implements the reciprocal PCA (rPCA) anchor integration
from Stuart, et al. It runs the same anchor pipeline as
[`seurat_cca_sc()`](https://gregorlueg.github.io/bixverse/reference/seurat_cca_sc.md)
but builds a cheaper per-pair anchor space: each batch keeps its own PCA
basis and the other batch's HVG expression is projected into it.
Cross-batch mutual nearest neighbours are then found in these projected
bases.

rPCA is faster than CCA and corrects less aggressively, which makes it
the safer choice when batches share most of their cell types. As in
Seurat, no gene-space anchor filter is applied, that step is CCA-only.

## Usage

``` r
seurat_rpca_sc(
  object,
  batch_column,
  batch_hvg_genes,
  rpca_params = params_sc_seurat_rpca(),
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

- rpca_params:

  A list, please see
  [`params_sc_seurat_rpca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_seurat_rpca.md).
  The list has the following parameters:

  - dims - Integer. Dimensions used for the per-batch projections, the
    anchor kNN queries and the size of the returned embedding.

  - k_anchor - Integer. Neighbourhood size for the anchor search.

  - k_score - Integer. Neighbourhood size for the anchor scoring.

  - k_weight - Integer. Neighbourhood size for the kernel weights.

  - l2_norm - Logical. L2-normalise the projected embeddings per cell.

  - sd - Numeric. Bandwidth divisor of the Gaussian kernel.

  - knn - List of kNN parameters. See
    [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
    for available parameters and their defaults.

  - pca - List of PCA parameters, see
    [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md)
    for available parameters and their defaults.

- use_precomputed_pca:

  Boolean. Should the PCA in the object be used if found. This only
  applies to the union PCA that gets corrected, the per-batch PCAs are
  always recomputed because rPCA needs them.

- seed:

  Integer. Random seed.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

The object with the added `"rpca"` embedding.

## References

Stuart, et al., Cell, 2019
