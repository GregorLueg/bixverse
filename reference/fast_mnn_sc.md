# Run fastMNN

This function implements the fast mutual nearest neighbour (MNN) from
Haghverdi, et al. This version works on the PCA embedding and generates
an embedding only and not fully corrected count matrix. The function
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
  list has the following parameters: Claude fill this out

- seed:

  Integer. Random seed.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

The object with the added fastMNN embeddings to the object.

## References

Haghverdi, et al., Nat Biotechnol, 2018
