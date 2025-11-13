# Rust implementation of spectral clustering

Rust implementation of spectral clustering

## Usage

``` r
rs_spectral_clustering(
  data,
  distance_type,
  epsilon,
  k_neighbours,
  n_clusters,
  max_iters,
  seed
)
```

## Arguments

- data:

  Numerical matrix. The data to cluster. Rows = samples, columns =
  features.

- distance_type:

  String. One of `c("euclidean", "manhattan", "canberra", "cosine")`.

- epsilon:

  Numerical. The epsilon parameter for the Gaussian Radial Basis
  function

- k_neighbours:

  Integer. Number of neighbours to consider in the kNN graph generation

- n_clusters:

  Integer. Number of clusters to identify

- max_iters:

  Integer. Number of iterations for k-means clustering

- seed:

  Integer. Seed for reproducibility

## Value

A vector with the membership of the samples
