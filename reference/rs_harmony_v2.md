# Harmony batch correction in Rust (version 2)

**\[experimental\]** This function implements the version 2 Harmony
algorithm from Patikas, et al., 2026.

## Usage

``` r
rs_harmony_v2(pca, harmony_params, batch_labels, seed, verbose)
```

## Arguments

- pca:

  Numerical matrix, i.e., the PCA matrix you want to correct.

- harmony_params:

  List. The parameters for the Harmony (v2) algorithm.

- batch_labels:

  List. Each element in the list needs to be a 0-indexed integer that
  represents the batch effects you wish to regress out.

- seed:

  Integer. Seed for reproducibility purposes.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

The batch-corrected Harmony (v2) embedding space.
