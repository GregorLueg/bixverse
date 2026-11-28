# Default parameters for the SCENIC ExtraTrees regression learner

Default parameters for the SCENIC ExtraTrees regression learner

## Usage

``` r
params_scenic_extra_trees_defaults()
```

## Value

A named list with the following elements:

- n_trees - Any. Number of trees to build. Defaults to `500L`.

- min_samples_leaf - Any. Minimum number of samples per leaf node.
  Defaults to `50L`.

- n_features_split - Any. Number of features considered per split. `0L`
  means `sqrt(n_features)`. Defaults to `0L`.

- n_thresholds - Any. Number of random thresholds tested per feature per
  node. Defaults to `1L`.

- max_depth - Any. Maximum tree depth. Defaults to `8L`.

- subsample_frac - Any. Optional. Fraction of cells to subsample per
  tree. `NULL` uses all cells. Defaults to `NULL`.
