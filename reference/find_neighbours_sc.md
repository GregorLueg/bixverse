# Find the neighbours for single cell.

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
  `hnsw`. Set `extract_knn = TRUE` in the kNN parameters to hand back
  the descent graph directly instead of beam searching it. That drops
  the query pass altogether, so it is much faster, but recall goes down
  a little.

- `ivf` - Inverted file index. Uses first k-means clustering to identify
  Voronoi cells and leverages these during querying. Works well on large
  data sets with high dimensionality and when you need to return large
  number of neighbours.

- `annoy` - Approximate nearest neighbours Oh Yeah. Tree-based index,
  used across different R single cell packages (Seurat, SCE). This
  version is purely memory-based.

Subsequently, the kNN is symmetrised into the neighbour graph the
clustering methods run on: the shared nearest neighbour graph (default),
the UMAP fuzzy simplicial set or the kNN union, see
[`params_knn_graph_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_graph_defaults.md).

## Usage

``` r
find_neighbours_sc(
  object,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  modality = c("rna", "adt"),
  neighbours_params = params_sc_neighbours(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` (or potentially other) class.

- embd_to_use:

  String. The embedding to use. Whichever you chose, it needs to be part
  of the object.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

- modality:

  String. One of `c("rna", "adt")`. You can only use `"adt"` on
  `SingleCellsMultiModal` class.

- neighbours_params:

  List. Output of
  [`params_sc_neighbours()`](https://gregorlueg.github.io/bixverse/reference/params_sc_neighbours.md).
  Carries the graph parameters (see
  [`params_knn_graph_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_graph_defaults.md))
  and the kNN parameters (see
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)).

- seed:

  Integer. For reproducibility.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

The object with the added kNN and neighbour graph.

## Examples

``` r
# kNN and the neighbour graph on top of the PCA
sc <- demo_single_cells(prepped = FALSE)
sc <- find_hvg_sc(sc, hvg_no = 30L, .verbose = FALSE)
sc <- calculate_pca_sc(sc, no_pcs = 10L, .verbose = FALSE)
sc <- find_neighbours_sc(
  sc,
  neighbours_params = params_sc_neighbours(knn = list(k = 15L)),
  .verbose = FALSE
)
dim(get_knn_mat(sc))
#> [1] 500  15

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
