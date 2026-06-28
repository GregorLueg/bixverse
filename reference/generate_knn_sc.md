# Generate a `SingleCellNearestNeighbour` from a single cell class

This function will generate the kNNs based on a given embedding.
Available algorithms are:

- `kmknn` - An exact kNN search that leverages k-means clustering under
  the hood to prune out data points. The default setting.

- `exhaustive` - An exhaustive, flat index. On smaller data sets often
  faster than the approximate nearest neighbour search algorithms.

- `hnsw` - Hierarchical Navigable Small World. A graph-based approximate
  nearest neighbour search algorithm; works well on large data sets. A
  benign race condition is leveraged during index build, making the
  build non-deterministic. Bigger impact on smaller data sets.

- `nndescent` - Nearest neighbour descent. Leverages concepts from
  `PyNNDescent` and works well on very large data sets similar to
  `hnsw`.

- `ivf` - Inverted file index. Uses first k-means clustering to identify
  Voronoi cells and leverages these during querying. Works well on large
  data sets with high dimensionality and when you need to return large
  number of neighbours.

- `annoy` - Approximate nearest neighbours Oh Yeah. Tree-based index,
  used across different R single cell packages (Seurat, SCE). This
  version is purely memory-based.

## Usage

``` r
generate_knn_sc(
  object,
  embd_to_use = "pca",
  cells_to_use = NULL,
  no_embd_to_use = NULL,
  neighbours_params = params_sc_neighbours(),
  seed = 42L,
  .validate_index = TRUE,
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
    only neighbours. Not used for this function.

  - pruning - Numeric. Weights below this threshold will be set to 0 in
    the generation of the sNN graph. Not used for this function.

  - snn_similarity - String. One of `c("rank", "jaccard")`. Defines how
    the weight from the SNN graph is calculated. For details, please see
    [`params_sc_neighbours()`](https://gregorlueg.github.io/bixverse/reference/params_sc_neighbours.md).
    Not used for this function.

  - knn - List of kNN parameters. See
    [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
    for available parameters and their defaults.

- seed:

  Integer. For reproducibility.

- .validate_index:

  Boolean. Shall an exhaustive search against a subset of cells be run
  to validate the approximate nearest neighbour index.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

Initialised `sc_knn` with the kNN data.
