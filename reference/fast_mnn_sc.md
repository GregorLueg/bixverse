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
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- batch_column:

  String. The column with the batch information in the obs data of the
  class.

- batch_hvg_genes:

  Integer vector. These are the highly variable genes, identified by a
  batch-aware method. Please refer to
  [`find_hvg_batch_aware_sc()`](find_hvg_batch_aware_sc.md) for more
  details.

- fastmnn_params:

  A list, please see [`params_sc_fastmnn()`](params_sc_fastmnn.md). The
  list has the following parameters:

  - sigma - Numeric. Bandwidth of the Gaussian smoothing kernel (as
    proportion of space radius).

  - cos_norm - Logical. Apply cosine normalisation before computing
    distances.

  - var_adj - Logical. Apply variance adjustment to avoid kissing
    effects.

  - no_pcs - Integer. Number of PCs to use for MNN calculations.

  - random_svd - Logical. Use randomised SVD.

  - knn - List of kNN parameters. See
    [`params_knn_defaults()`](params_knn_defaults.md) for available
    parameters and their defaults.

- seed:

  Integer. Random seed.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

The object with the added fastMNN embeddings to the object.

## References

Haghverdi, et al., Nat Biotechnol, 2018
