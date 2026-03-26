# Generates the kNN graph with additional distances

This function is a wrapper over the Rust-based generation of the
approximate nearest neighbours.

## Usage

``` r
rs_sc_knn_w_dist(embd, knn_params, verbose, seed)
```

## Arguments

- embd:

  Numerical matrix. The embedding matrix to use to generate the kNN
  graph.

- knn_params:

  List. The kNN parameters defined by
  [`params_sc_neighbours()`](https://gregorlueg.github.io/bixverse/reference/params_sc_neighbours.md).

- verbose:

  Boolean. Controls verbosity of the function and returns how long
  certain operations took.

- seed:

  Integer. Seed for reproducibility purposes.

## Value

A list with:

- indices - An integer matrix representing the indices of the
  approximate nearest neighbours.

- dist - An numerical matrix representing the distances to the nearest
  neighbours.

- dist_metric - String representing the used distance metric.
