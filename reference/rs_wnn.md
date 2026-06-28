# Run the weighted nearest neighbour algorithm

**\[experimental\]** This provides a Rust-based implementation of the
WNN algorithm from Hao, et al.

## Usage

``` r
rs_wnn(modality_emb_one, modality_emb_two, wnn_params, seed, verbose)
```

## Arguments

- modality_emb_one:

  Numerical matrix of the first modality. For example the PCA (or other
  embeddings) from the transcriptomics.

- modality_emb_two:

  Numerical matrix of the second modality. For example the PCA (or other
  embeddings) from the ADT counts.

- wnn_params:

  Named list. The weighted nearest neighbour parameters.

- seed:

  Integer. For reproducibility purposes.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with

- indices An integer matrix representing the indices of the approximate
  nearest neighbours.

- dist - An numerical matrix representing the distances to the nearest
  neighbours.

- modality_one_weights - The weights of the first modality.

- modality_two_weights - The weights of the second modality.

## References

Hao et al., Cell, 2021
