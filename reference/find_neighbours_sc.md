# Find the neighbours for single cell.

This function will generate the kNNs based on a given embedding. Three
different algorithms are implemented with different speed and accuracy
to approximate the nearest neighbours. `"annoy"` is more rapid and based
on the `Approximate Nearest Neigbours Oh Yeah` algorithm; `"hnsw"`
implements a `Hierarchical Navigatable Small Worlds` vector search that
is slower, but more precise. Lastly, there is the option of
`"nndescent"`, a Rust-based implementation of the PyNNDescent algorithm.
This version skips the index generation and can be faster on smaller
data sets. Subsequently, the kNN data will be used to generate an sNN
igraph for clustering methods.

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

  `single_cell_exp` class.

- embd_to_use:

  String. The embedding to use. Whichever you chose, it needs to be part
  of the object.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

- neighbours_params:

  List. Output of [`params_sc_neighbours()`](params_sc_neighbours.md). A
  list with the following items:

  - k - Integer. Number of neighbours to identify.

  - knn_algorithm - String. One of `c("annoy", "hnsw", "nndescent")`.
    `"hnsw"` takes longer, is more precise and more memory friendly.
    `"annoy"` is faster, less precise and will take more memory.
    `"nndescent"` skips index generation and can be faster on small
    datasets.

  - n_trees - Integer. Number of trees to use for the `annoy` algorithm.
    The higher, the longer the algorithm takes, but the more precise the
    approximated nearest neighbours.

  - search_budget - Integer. Search budget per tree for the `annoy`
    algorithm. The higher, the longer the algorithm takes, but the more
    precise the approximated nearest neighbours.

  - ann_dist - String. One of `c("cosine", "euclidean")`.

  - max_iter - Integer. Maximum iterations for the `"nndescent"` method.

  - rho - Numeric. Sampling rate for the `"nndescent"` method.

  - delta - Numeric. Early termination criterium for the `"nndescent"`
    method.

  - full_snn - Boolean. Shall the sNN graph be generated across all
    cells (standard in the `bluster` package.) Defaults to `FALSE`.

  - pruning - Value below which the weight in the sNN graph is set to 0.

  - snn_similarity - String. One of `c("rank", "jaccard")`. Defines how
    the weight form the SNN graph is calculated. For details, please see
    [`params_sc_neighbours()`](params_sc_neighbours.md).

- seed:

  Integer. For reproducibility.

- .verbose:

  Boolean. Controls verbosity and returns run times.

## Value

The object with added KNN matrix.
