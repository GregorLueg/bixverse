# Default parameters for the SCENIC RandomForest regression learner

Default parameters for the SCENIC RandomForest regression learner

## Usage

``` r
params_scenic_random_forest_defaults()
```

## Value

A named list with the following elements:

- n_trees - Any. Number of trees to build. Defaults to `250L`.

- min_samples_leaf - Any. Minimum number of samples per leaf node.
  Defaults to `50L`.

- n_features_split - Any. Number of features considered per split. `0L`
  means `sqrt(n_features)`. Defaults to `0L`.

- subsample_rate - Any. Fraction of samples drawn per tree, without
  replacement unless `bootstrap = TRUE`. Defaults to `0.632`.

- bootstrap - Any. Whether to sample with replacement. Defaults to
  `FALSE`.

- max_depth - Any. Maximum tree depth. Defaults to `8L`.

- subsample_frac - Any. Optional. Fraction of cells to subsample per
  tree. If set, overrides `subsample_rate`. Defaults to `NULL`.
