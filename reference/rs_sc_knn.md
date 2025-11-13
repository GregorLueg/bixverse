# Generates the kNN graph

This function is a wrapper over the Rust-based generation of the
approximate nearest neighbours. You have several options to get the
approximate nearest neighbours:

- `"annoy"`: leverages binary trees to generate rapidly in a parallel
  manner an index. Good compromise of index generation, querying speed.

- `"hnsw"`: uses a hierarchical navigatable small worlds index under the
  hood. The index generation takes more long, but higher recall and
  ideal for very large datasets due to subdued memory pressure.

- `"nndescent"`: an index-free approximate nearest neighbour algorithm
  that is ideal for small, ephemeral kNN graphs.

## Usage

``` r
rs_sc_knn(embd, knn_params, verbose, seed)
```

## Arguments

- embd:

  Numerical matrix. The embedding matrix to use to generate the kNN
  graph.

- knn_params:

  List. The kNN parameters defined by
  [`params_sc_neighbours()`](params_sc_neighbours.md).

- verbose:

  Boolean. Controls verbosity of the function and returns how long
  certain operations took.

- seed:

  Integer. Seed for reproducibility purposes.

## Value

A integer matrix of N x k with N being the number of cells and k the
number of neighbours.
