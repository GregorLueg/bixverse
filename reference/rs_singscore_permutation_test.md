# Rust version of the singscore permutation test

**\[experimental\]** For `n_permutations` iterations, draws random gene
indices of the same total size as the real gene set(s), scores them with
the same options as the real call, and builds a per-sample null
distribution. Returns empirical one-tailed p-values:
`max(1 / n_permutations, mean(null > observed))`.

## Usage

``` r
rs_singscore_permutation_test(
  ranks,
  up_set,
  down_set,
  center_score,
  known_direction,
  stable,
  n_permutations,
  seed
)
```

## Arguments

- ranks:

  Numerical matrix. The ranked expression matrix.

- up_set:

  Integer vector. Zero-indexed.

- down_set:

  Integer vector or NULL.

- center_score, known_direction, stable:

  Booleans. Should match the values used for the real
  [`rs_singscore_single()`](https://gregorlueg.github.io/bixverse/reference/rs_singscore_single.md)
  call.

- n_permutations:

  Integer. Number of random draws (B).

- seed:

  Integer. RNG seed.

## Value

A named list with `observed_scores` (length n_samples),
`null_distribution` (B × n_samples matrix), and `p_values` (length
n_samples).
