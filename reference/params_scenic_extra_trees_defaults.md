# Default parameters for the SCENIC ExtraTrees regression learner

Default parameters for the SCENIC ExtraTrees regression learner

## Usage

``` r
params_scenic_extra_trees_defaults()
```

## Value

A list with the following parameters:

- n_trees - Integer. Number of trees to build. Defaults to `500L`.

- min_samples_leaf - Integer. Minimum number of samples required at a
  leaf node. Defaults to `50L`.

- n_features_split - Integer. Number of features considered per split.
  `0L` resolves to `sqrt(n_features)` at runtime. Defaults to `0L`.

- n_thresholds - Integer. Number of random thresholds to evaluate per
  feature per node. Defaults to `1L`.

- max_depth - Integer. Maximum depth of each tree. Defaults to `8L`.

- subsample_frac - Optional numeric. Fraction of cells to subsample per
  tree. Defaults to `NULL`.
