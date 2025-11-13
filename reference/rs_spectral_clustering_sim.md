# Rust implementation of spectral clustering

Spectral clustering on a pre-calculated similarity matrix.

## Usage

``` r
rs_spectral_clustering_sim(
  similarities,
  k_neighbours,
  n_clusters,
  max_iters,
  seed
)
```

## Arguments

- similarities:

  Numerical matrix representing the similarities. Needs to be symmetric!

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
