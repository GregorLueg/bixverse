# SCENIC on MetaCells

Assumes that the sparse data is pre-filtered for the genes you wish to
include. The indices need to be 0-indexed.

## Usage

``` r
rs_mc_scenic(sparse_data, tf_indices, scenic_params, seed, verbose)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `format`.

- tf_indices:

  Integer vector. The indices of the transcription factors.

- scenic_params:

  Named list. Contains all of the parameters need for SCENIC.

- seed:

  Integer. Controls reproducibility of the function.

- verbose:

  Boolean. Controls the verbosity of the function.

## Value

A gene x TF importance matrix
