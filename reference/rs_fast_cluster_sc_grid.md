# Runs fast Louvain cluster on the data (with multiple seeds)

**\[experimental\]** Runs first k-means clustering, followed by a kNN
detection on the centroids to then run Louvain clustering with several
seeds (based on the original one) on the graph and propagate the
membership back to the original data. Returns additional metrics around
cluster stability and community conductance.

## Usage

``` r
rs_fast_cluster_sc_grid(
  embd,
  km_type,
  resolutions,
  n_centroids,
  fc_params,
  snn,
  return_kmeans,
  no_seeds,
  seed,
  verbose
)
```

## Arguments

- embd:

  Numeric matrix. The original embedding.

- km_type:

  String. One of `c("kmeans", "minibatch")` for the type of k means
  clustering to run.

- resolutions:

  Numeric vector. The Louvain resolutions to iterate through.

- n_centroids:

  Optional integer. The number of k-means centroids. If not provided,
  defaults to `floor(sqrt(nrow(embd)))`.

- fc_params:

  Named list. The fast clustering parameters.

- snn:

  Boolean. Shall the kNN graph be additionally transformed into an sNN
  graph.

- return_kmeans:

  Boolean. Shall the k-means centroid assignments be returned alongside
  the memberships.

- no_seeds:

  Integer. Number of additional seeds to use. Should be \>=2.

- seed:

  Integer. For reproducibility.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following elements:

- membership - A list with `memberships` (one integer vector per
  resolution, from the seed with the best conductance) and `stats` (list
  with `mean_ari`, `median_ari`, `mean_conductance`,
  `median_conductance` and `mean_n_comms`, one value per resolution).

- k_means_cluster - Integer vector with the k-means cluster per cell if
  `return_kmeans = TRUE`, otherwise `NULL`.

- centroids - Numerical matrix with the k-means centroids if
  `return_kmeans = TRUE`, otherwise `NULL`.
