# Run ForceAtlas2 on a SingleCells/MetaCells object

Wrapper around
[`manifoldsR::forceatlas2()`](https://gregorlueg.github.io/manifoldsR/reference/forceatlas2.html)
and
[`manifoldsR::forceatlas2_from_graph()`](https://gregorlueg.github.io/manifoldsR/reference/forceatlas2_from_graph.html)
for the `SingleCells` and `MetaCells` classes. ForceAtlas2 is a
force-directed graph layout: every edge pulls, every pair of cells
pushes, and gravity keeps disconnected components from drifting off. It
is what scanpy's `draw_graph` does.

`graph` picks what gets laid out:

- `"knn"` - The kNN graph (cached if `use_knn = TRUE`, otherwise built
  from the chosen embedding), turned into the UMAP fuzzy union graph
  first. This matches scanpy.

- `"snn"` - The sNN graph from
  [`find_neighbours_sc()`](https://gregorlueg.github.io/bixverse/reference/find_neighbours_sc.md),
  i.e. the same graph the Leiden/Louvain clustering runs on. `use_knn`,
  `embd_to_use`, `no_embd_to_use`, `k`, `knn_method` and `nn_params` do
  not apply here, and neither do the graph and initialisation knobs in
  `fa2_params`. Pass `init_embd` to start from an existing 2D embedding,
  otherwise the layout starts from random positions.

## Usage

``` r
forceatlas2_sc(
  object,
  graph = c("knn", "snn"),
  use_knn = TRUE,
  embd_to_use = "pca",
  slot_name = "fa2",
  no_embd_to_use = NULL,
  init_embd = NULL,
  modality = c("rna", "adt", "wnn"),
  k = 15L,
  knn_method = c("kmknn", "hnsw", "balltree", "annoy", "nndescent", "exhaustive"),
  nn_params = manifoldsR::params_nn(),
  fa2_params = manifoldsR::params_fa2(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` class.

- graph:

  String. Which graph to lay out. One of `c("knn", "snn")`. Defaults to
  `"knn"`.

- use_knn:

  Boolean. Use the kNN graph found in the object. Defaults to `TRUE`. If
  not available, will default to the embedding. `"knn"` only.

- embd_to_use:

  String. The embedding to build the kNN graph from. Must be available
  in the object. `"knn"` only.

- slot_name:

  String. The name of this embedding within the object. Defaults to
  `"fa2"`.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used. `"knn"` only.

- init_embd:

  Optional string. Name of a stored 2D embedding, e.g. `"umap"`, to
  initialise the layout with. `"snn"` only. Defaults to `NULL`.

- modality:

  String. On which modality to run ForceAtlas2. One of
  `c("rna", "adt", "wnn")`. The two latter options are only available
  for multi-modal versions with the added data.

- k:

  Integer. Number of nearest neighbours. Defaults to `15L`. `"knn"`
  only.

- knn_method:

  String. Approximate nearest neighbour algorithm. One of `"hnsw"`,
  `"balltree"`, `"annoy"`, `"nndescent"`, or `"exhaustive"`. `"knn"`
  only.

- nn_params:

  Named list. See
  [`manifoldsR::params_nn()`](https://gregorlueg.github.io/manifoldsR/reference/params_nn.html).
  `"knn"` only.

- fa2_params:

  Named list. See
  [`manifoldsR::params_fa2()`](https://gregorlueg.github.io/manifoldsR/reference/params_fa2.html).

- seed:

  Integer. For reproducibility.

- .verbose:

  Boolean. Controls verbosity.

## Value

The object with a `"fa2"` embedding added.

## References

Jacomy, et al., PLoS ONE, 2014

## Examples

``` r
# ForceAtlas2 on the cached kNN graph, then on the sNN graph
sc <- demo_single_cells()
sc <- forceatlas2_sc(sc, .verbose = FALSE)
sc <- forceatlas2_sc(
  sc,
  graph = "snn",
  slot_name = "fa2_snn",
  .verbose = FALSE
)
dim(get_embedding(sc, "fa2_snn"))
#> [1] 500   2

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
