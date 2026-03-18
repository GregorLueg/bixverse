# Default parameters for the SCENIC GradientBoosting (GRNBoost2) regression learner

Default parameters for the SCENIC GradientBoosting (GRNBoost2)
regression learner

## Usage

``` r
params_scenic_gradient_boosting_defaults()
```

## Value

A list with the following parameters:

- n_trees_max - Integer. Maximum number of boosting rounds. Early
  stopping usually triggers well before this limit. Defaults to `1000L`.

- learning_rate - Numeric. Shrinkage applied to each tree's predictions.
  Defaults to `0.01`.

- max_depth - Integer. Maximum depth of each tree. Shallow trees (3-5)
  work best for GBM. Defaults to `3L`.

- min_samples_leaf - Integer. Minimum number of training samples
  required at a leaf node. Defaults to `50L`.

- early_stop_window - Integer. Number of recent OOB improvements to
  average for the early stopping criterion. Stops when the rolling
  average drops to zero or below. Defaults to `25L`.

- subsample_rate - Numeric. Fraction of samples used for training each
  tree. The complement forms the OOB set. Defaults to `0.9`.

- n_features_split - Integer. Number of features to evaluate per split.
  `0L` means all features (recommended with histogram subtraction).
  Defaults to `0L`.
