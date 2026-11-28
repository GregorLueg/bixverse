# Default parameters for the SCENIC GradientBoosting (GRNBoost2) regression learner

Default parameters for the SCENIC GradientBoosting (GRNBoost2)
regression learner

## Usage

``` r
params_scenic_gradient_boosting_defaults()
```

## Value

A named list with the following elements:

- n_trees_max - Any. Maximum number of boosting rounds. Early stopping
  usually triggers well before this. Defaults to `1000L`.

- learning_rate - Any. Shrinkage applied to each tree's predictions.
  Defaults to `0.01`.

- max_depth - Any. Maximum tree depth. Shallow trees (3 to 5) work best.
  Defaults to `3L`.

- min_samples_leaf - Any. Minimum number of training samples per leaf
  node. Defaults to `50L`.

- early_stop_window - Any. Number of recent out-of-bag improvements
  averaged for early stopping. Boosting stops once that rolling average
  drops to zero or below. Defaults to `25L`.

- subsample_rate - Any. Fraction of samples used to train each tree. The
  rest form the out-of-bag set for early stopping. Defaults to `0.9`.

- n_features_split - Any. Number of features evaluated per split. `0L`
  means all features (unlike the forest learners, where it means
  `sqrt(n_features)`). Defaults to `0L`.
