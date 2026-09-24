# Runs fast Louvain cluster on the data

**\[experimental\]** Runs first k-means clustering, followed by a kNN
detection on the centroids to then run Louvain clustering on the graph
and propagate the membership back to the original data.

## Usage

``` r
rs_fast_cluster_sc(
  embd,
  km_type,
  resolutions,
  n_centroids,
  fc_params,
  snn,
  return_kmeans,
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

- seed:

  Integer. For reproducibility.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following elements:

- membership - List with one integer membership vector per resolution.

- k_means_cluster - Integer vector with the k-means cluster per cell if
  `return_kmeans = TRUE`, otherwise `NULL`.

- centroids - Numerical matrix with the k-means centroids if
  `return_kmeans = TRUE`, otherwise `NULL`.
