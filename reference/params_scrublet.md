# Wrapper function for Scrublet doublet detection parameters

Constructor for the various Scrublet parameters. In this case, the
default for the kNN graph generation was set to `"hnsw"` as this
algorithm showed the best performance in different empirical benchmarks.

## Usage

``` r
params_scrublet(
  sim_doublet_ratio = 1.5,
  expected_doublet_rate = 0.1,
  stdev_doublet_rate = 0.02,
  n_bins = 100L,
  manual_threshold = NULL,
  normalisation = list(),
  hvg = list(),
  pca = list(),
  knn = list(k = 0L, knn_method = "hnsw")
)
```

## Arguments

- sim_doublet_ratio:

  Numeric. Number of doublets to simulate relative to the number of
  observed cells. For example, 2.0 simulates twice as many doublets as
  there are cells. Defaults to `1.5`.

- expected_doublet_rate:

  Numeric. Expected doublet rate for the experiment, typically 0.05-0.10
  depending on cell loading. Must be between 0 and 1. Defaults to `0.1`.

- stdev_doublet_rate:

  Numeric. Uncertainty in the expected doublet rate. Defaults to `0.02`.

- n_bins:

  Integer. Number of bins for histogram-based automatic threshold
  detection. Typically 50-100. Defaults to `100L`.

- manual_threshold:

  Optional numeric. Manual doublet score threshold. If `NULL` (default),
  threshold is automatically detected from simulated doublet score
  distribution.

- normalisation:

  List. Optional overrides for normalisation parameters. See
  [`params_norm_doublet_detection_defaults()`](params_norm_doublet_detection_defaults.md)
  for available parameters: `log_transform`, `mean_center`,
  `normalise_variance`, `target_size`.

- hvg:

  List. Optional overrides for highly variable gene selection
  parameters. See [`params_hvg_defaults()`](params_hvg_defaults.md) for
  available parameters: `min_gene_var_pctl`, `hvg_method`, `loess_span`,
  `clip_max`.

- pca:

  List. Optional overrides for PCA parameters. See
  [`params_pca_defaults()`](params_pca_defaults.md) for available
  parameters: `no_pcs`, `random_svd`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](params_knn_defaults.md) for available
  parameters: `k`, `knn_method`, `dist_metric`, `search_budget`,
  `n_trees`, `nn_max_iter`, `rho`, `delta`.

## Value

A named list with all Scrublet parameters, combining defaults with any
user-specified overrides.
