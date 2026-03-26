# SCENIC: Select TF-gene pairs by per-gene importance threshold

For each gene (row), computes mean + n_sd \* SD of the importance scores
across all TFs and retains only pairs exceeding that threshold.

## Usage

``` r
rs_importance_threshold(matrix, n_sd, min_value)
```

## Arguments

- matrix:

  Numeric matrix with genes (rows) x TFs (columns) importance values.

- n_sd:

  Float. Number of standard deviations above the mean to use as the
  per-gene threshold.

- min_value:

  Optional float. Absolute minimum importance score. Pairs below this
  are excluded even if they pass the per-gene threshold.

## Value

A list with three vectors: tf, gene, importance
