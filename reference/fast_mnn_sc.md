# Run fastMNN

This function implements the fast mutual nearest neighbour (MNN) from
Haghverdi, et al. This version works on the PCA embedding and generates
an embedding only and not a fully corrected count matrix. The function
will iterate through the batches, identify the MNN and generate
correction vectors and generate a corrected embedding which is added to
the function.

## Usage

``` r
fast_mnn_sc(
  object,
  batch_column,
  batch_hvg_genes,
  fastmnn_params = params_sc_fastmnn(),
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

- fastmnn_params:

  A list, please see
  [`params_sc_fastmnn()`](https://gregorlueg.github.io/bixverse/reference/params_sc_fastmnn.md).
  The list has the following parameters:

  - sigma - Numeric. Bandwidth of the Gaussian smoothing kernel (as
    proportion of space radius).

  - cos_norm - Logical. Apply cosine normalisation before computing
    distances.

  - var_adj - Logical. Apply variance adjustment to avoid kissing
    effects.

  - no_pcs - Integer. Number of PCs to use for MNN calculations.

  - random_svd - Logical. Use randomised SVD.

  - knn - List of kNN parameters. See
    [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
    for available parameters and their defaults.

- use_precomputed_pca:

  Boolean. Should the PCA in the object be used if found. If you decide
  to do this, make sure that you have run the PCA on the batch-aware HVG
  ideally.

- seed:

  Integer. Random seed.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

The object with the added fastMNN embeddings to the object.

## References

Haghverdi, et al., Nat Biotechnol, 2018
