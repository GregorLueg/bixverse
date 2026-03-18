# Default parameters for the SCENIC RandomForest regression learner

Default parameters for the SCENIC RandomForest regression learner

## Usage

``` r
params_scenic_random_forest_defaults()
```

## Value

A list with the following parameters:

- n_trees - Integer. Number of trees to build. Defaults to `250L`.

- min_samples_leaf - Integer. Minimum number of samples required at a
  leaf node. Defaults to `50L`.

- n_features_split - Integer. Number of features considered per split.
  `0L` resolves to `sqrt(n_features)` at runtime. Defaults to `0L`.

- subsample_rate - Numeric. Fraction of samples to draw per tree.
  Defaults to `0.632`.

- bootstrap - Logical. Whether to sample with replacement. Defaults to
  `FALSE`.

- max_depth - Integer. Maximum depth of each tree. Defaults to `8L`.

- subsample_frac - Optional numeric. Fraction of cells to subsample per
  tree. If set, overrides `subsample_rate`. Defaults to `NULL`.
