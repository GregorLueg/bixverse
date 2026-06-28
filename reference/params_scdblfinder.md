# Wrapper function for scDblFinder doublet detection parameters

Constructor for the scDblFinder parameters. This method combines
cluster-aware doublet simulation with a gradient-boosted classifier
trained on engineered features.

## Usage

``` r
params_scdblfinder(
  n_genes = 1352L,
  doublet_ratio = 1,
  heterotypic_bias = 1,
  cluster_resolution = 1,
  cluster_iters = 10L,
  fast_cluster = FALSE,
  n_iterations = 3L,
  n_trees = 200L,
  max_depth = 4L,
  learning_rate = 0.3,
  min_samples_leaf = 20L,
  subsample_rate = 0.75,
  cv_folds = 5L,
  cv_early_stop = 2L,
  se_fraction = 1,
  include_pcs = 19L,
  expected_doublet_rate = NULL,
  cxds_genes = NULL,
  manual_threshold = NULL,
  normalisation = list(mean_center = TRUE),
  pca = list(),
  knn = list(k = 0L),
  fast_cluster_params = list()
)
```

## Arguments

- n_genes:

  Integer. Number of top-expressed genes to use as features. Defaults to
  `1352L`.

- doublet_ratio:

  Numeric. Ratio of simulated doublets to observed cells. Defaults to
  `1.0`.

- heterotypic_bias:

  Numeric. Fraction of simulated pairs forced to come from different
  clusters (0-1). Defaults to `1.0`.

- cluster_resolution:

  Numeric. Resolution for the initial Louvain clustering. Defaults to
  `1.0`.

- cluster_iters:

  Integer. Number of Louvain iterations per clustering step. Defaults to
  `10L`.

- fast_cluster:

  Boolean. Shall fast Louvain clustering be applied, i.e., k-means
  clustering and use the centroids for kNN graph generation and Louvain
  clustering with then backpropagating the membership based on centroid
  proximity.

- n_iterations:

  Integer. Number of refinement iterations. Typically 2-3. Defaults to
  `3L`.

- n_trees:

  Integer. Maximum number of boosting rounds for the GBM classifier.
  Defaults to `200L`.

- max_depth:

  Integer. Maximum tree depth. Shallow trees (3-5) work best. Defaults
  to `4L`.

- learning_rate:

  Numeric. Shrinkage applied to each tree. Defaults to `0.3`.

- min_samples_leaf:

  Integer. Minimum training samples per leaf. Defaults to `20L`.

- subsample_rate:

  Numeric. Fraction of samples used per tree. Defaults to `0.75`.

- cv_folds:

  Integer. Number of cross-validation folds for boosting round
  selection. Defaults to `5L`.

- cv_early_stop:

  Integer. Early stopping patience per CV fold. Defaults to `2L`.

- se_fraction:

  Numeric. Multiplier on the standard error for the SE rule used in
  round selection. Defaults to `1.0`

- include_pcs:

  Integer. Number of leading principal components to include as
  classifier features. Defaults to `19L`.

- expected_doublet_rate:

  Optional numeric. Expected doublet rate as a percentage. If not
  provided, will be calculated internally.

- cxds_genes:

  Optional integer. Number of CXDS genes to consider. If not provided,
  defaults to `500L`.

- manual_threshold:

  Optional numeric. Manual score threshold. If `NULL` (default),
  expected-rate thresholding is used.

- normalisation:

  List. Optional overrides for normalisation parameters. See
  [`params_norm_doublets_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_norm_doublets_defaults.md).

- pca:

  List. Optional overrides for PCA parameters. See
  [`params_pca_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_pca_defaults.md).

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md).
  NNDescent works better for the larger k-values often used here.

- fast_cluster_params:

  List. Optional overrides for the fast clustering parameters. Only
  relevant if `fast_cluster = TRUE`. See
  [`params_fast_cluster_default()`](https://gregorlueg.github.io/bixverse/reference/params_fast_cluster_default.md)
  for available parameters: `km_type`, `n_centroids`, `kmeans_iters` and
  `batch_size`.

## Value

A named list with all scDblFinder parameters.
