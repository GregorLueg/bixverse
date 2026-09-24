# Run the weighted nearest neighbour algorithm

**\[experimental\]** This provides a Rust-based implementation of the
WNN algorithm from Hao, et al. Both embeddings are L2-normalised per
cell, a kNN graph with `knn_range` neighbours is built per modality, and
the per-cell modality weights are then used to fuse both into one kNN
graph.

## Usage

``` r
rs_wnn(modality_emb_one, modality_emb_two, wnn_params, seed, verbose)
```

## Arguments

- modality_emb_one:

  Numerical matrix of the first modality, cells x dimensions. For
  example the PCA (or other embeddings) from the transcriptomics.

- modality_emb_two:

  Numerical matrix of the second modality, same cells in the same row
  order. For example the PCA (or other embeddings) from the ADT counts.

- wnn_params:

  Named list. The weighted nearest neighbour parameters.

- seed:

  Integer. For reproducibility purposes.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with

- indices - Integer matrix of cells x neighbours with the indices
  (0-indexed!) of the weighted nearest neighbours.

- dist - Numerical matrix with the distances to these neighbours.

- dist_metric - String. Always `"kernelised pseudo-distance"`.

- modality_one_weights - Numerical vector with the per-cell weights of
  the first modality.

- modality_two_weights - Numerical vector with the per-cell weights of
  the second modality.

## References

Hao et al., Cell, 2021
