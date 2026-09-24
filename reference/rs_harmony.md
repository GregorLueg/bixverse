# Harmony batch correction in Rust

**\[experimental\]** This function implements the Harmony algorithm from
Korsunsky et al., 2019.

## Usage

``` r
rs_harmony(pca, harmony_params, batch_labels, seed, verbose)
```

## Arguments

- pca:

  Numerical matrix, i.e., the PCA matrix you want to correct.

- harmony_params:

  List. The parameters for the Harmony algorithm.

- batch_labels:

  List. Each element needs to be a 0-indexed integer vector, one per
  batch variable you wish to regress out.

- seed:

  Integer. Seed for reproducibility purposes.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

Numerical matrix, cells x dimensions, with the batch-corrected Harmony
embedding.
