# Generates the kNN graph

**\[experimental\]** This function is a wrapper over the Rust-based
generation of the approximate nearest neighbours.

## Usage

``` r
rs_sc_knn(embd, knn_params, validate_index, verbose, seed)
```

## Arguments

- embd:

  Numerical matrix. The embedding matrix to use to generate the kNN
  graph.

- knn_params:

  List. The kNN parameters defined by
  [`params_sc_neighbours()`](https://gregorlueg.github.io/bixverse/reference/params_sc_neighbours.md).

- validate_index:

  Boolean. If you want to validate the index via an exhaustive search in
  a subset of cells.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

- seed:

  Integer. Seed for reproducibility purposes.

## Value

A integer matrix of N x k with N being the number of cells and k the
number of neighbours.
