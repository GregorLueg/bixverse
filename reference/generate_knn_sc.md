# Generate the KNN data with distances

This function will generate the kNNs based on a given embedding. Three
different algorithms are implemented with different speed and accuracy
to approximate the nearest neighbours. `"annoy"` is more rapid and based
on the `Approximate Nearest Neigbours Oh Yeah` algorithm; `"hnsw"`
implements a `Hierarchical Navigatable Small Worlds` vector search that
is slower, but more precise. Lastly, there is the option of
`"nndescent"`, a Rust-based implementation of the PyNNDescent algorithm.
This version skips the index generation and can be faster on smaller
data sets. This version of the function returns an `sc_knn` object that
can be used in other functions.

## Usage

``` r
generate_knn_sc(
  object,
  embd_to_use = "pca",
  cells_to_use = NULL,
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

- cells_to_use:

  String. Optional cell names to include in the generation of the kNN
  graph. If `NULL` all (filtered) cells in the object will be used.

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

Initialised `sc_knn` with the kNN data.
