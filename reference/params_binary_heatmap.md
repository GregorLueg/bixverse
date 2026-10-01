# Wrapper function for parameters for binary heatmap data

Controls filtering, ordering and binning in
[`extract_binary_heatmap_data()`](https://gregorlueg.github.io/bixverse/reference/extract_binary_heatmap_data.md).
Features are clustered on the Jaccard distance, samples on the Hamming
distance within their group. Groups larger than `max_cluster_n` are
ordered by barycentre instead, which skips the distance matrix. More
samples than `max_cols` get binned within their group into fraction-on
columns.

## Usage

``` r
params_binary_heatmap(
  min_frac_on = 0.01,
  max_frac_on = 0.99,
  cluster_features = TRUE,
  cluster_samples = TRUE,
  max_cols = 2000L,
  max_cluster_n = 2000L
)
```

## Arguments

- min_frac_on:

  Numeric. Features on in a smaller fraction of samples are dropped.
  Defaults to `0.01`.

- max_frac_on:

  Numeric. Features on in a larger fraction of samples are dropped.
  Defaults to `0.99`.

- cluster_features:

  Boolean. Shall the features be clustered within their group. Defaults
  to `TRUE`.

- cluster_samples:

  Boolean. Shall the samples be ordered within their group. Defaults to
  `TRUE`.

- max_cols:

  Integer. Maximum number of columns before samples get binned. Defaults
  to `2000L`.

- max_cluster_n:

  Integer. Groups with more samples than this are ordered by barycentre
  instead of hierarchical clustering. Defaults to `2000L`.

## Value

A named list with the following elements:

- min_frac_on - Numeric. Features on in a smaller fraction of samples
  are dropped. Defaults to `0.01`.

- max_frac_on - Numeric. Features on in a larger fraction of samples are
  dropped. Defaults to `0.99`.

- cluster_features - Boolean. Shall the features be clustered within
  their group. Defaults to `TRUE`.

- cluster_samples - Boolean. Shall the samples be ordered within their
  group. Defaults to `TRUE`.

- max_cols - Integer. Maximum number of columns before samples get
  binned. Defaults to `2000L`.

- max_cluster_n - Integer. Groups with more samples than this are
  ordered by barycentre instead of hierarchical clustering. Defaults to
  `2000L`.
