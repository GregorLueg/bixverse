# Run MELD

**\[experimental\]** This implements a Rust-based version of the MELD
algorithm, see Burkhardt, et al. Nat. Biotechnol., 2021.

## Usage

``` r
rs_meld_sc(embd, knn_data, meld_params, labels, n_labels, seed, verbose)
```

## Arguments

- embd:

  Numeric matrix. The original embedding that was used to generate the
  kNN graph.

- knn_data:

  Optional named list. This contains pre-computed kNN data (including
  distances). The user has to ensure consistency! If provided, this will
  be used.

- meld_params:

  Named list. Contains the parameters to use for MELD.

- labels:

  Integer. The labels of the different groups. (1-indexed!)

- n_labels:

  Integer. Number of labels represented in the data.

- seed:

  Integer. For reproducibility.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- raw_scores - The raw MELD scores

- norm_scores - Negative values were clamped to 0 and the rows L1
  normalised. This yields probability-like values.
