# Find the neighbours for single cell.

This function will generate the kNNs based on a given embedding. Four
different algorithms are implemented with different speed and accuracy
to approximate the nearest neighbours. `"hnsw"` implements a
Hierarchical Navigatable Small Worlds vector search that has slower
index generation but high precision. `"annoy"` is based on the
Approximate Nearest Neighbours Oh Yeah algorithm and is more rapid in
terms of index generation, but querying on large data sets can be slow.
`"nndescent"` is a Rust-based implementation of the PyNNDescent
algorithm and is a good all-rounder and performs well on very large data
sets. `"exhaustive"` performs exact nearest neighbour search.
Subsequently, the kNN data will be used to generate an sNN igraph for
clustering methods.

## Usage

``` r
find_neighbours_sc(
  object,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  neighbours_params = params_sc_neighbours(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- embd_to_use:

  String. The embedding to use. Whichever you chose, it needs to be part
  of the object.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

- neighbours_params:

  List. Output of
  [`params_sc_neighbours()`](https://gregorlueg.github.io/bixverse/reference/params_sc_neighbours.md).
  A list with the following items:

  - full_snn - Boolean. Shall the full shared nearest neighbour graph be
    generated that generates edges between all cells instead of between
    only neighbours.

  - pruning - Numeric. Weights below this threshold will be set to 0 in
    the generation of the sNN graph.

  - snn_similarity - String. One of `c("rank", "jaccard")`. Defines how
    the weight from the SNN graph is calculated. For details, please see
    [`params_sc_neighbours()`](https://gregorlueg.github.io/bixverse/reference/params_sc_neighbours.md).

  - knn - List of kNN parameters. See
    [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
    for available parameters and their defaults.

- seed:

  Integer. For reproducibility.

- .verbose:

  Boolean. Controls verbosity and returns run times.

## Value

The object with added KNN matrix.
