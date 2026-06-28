# Run UMAP on a SingleCells/MetaCells object

Wrapper around
[`manifoldsR::umap()`](https://gregorlueg.github.io/manifoldsR/reference/umap.html)
for the `SingleCells` and `MetaCells` classes. UMAP produces a
low-dimensional embedding that emphasises local neighbourhood structure
while being computationally efficient via its negative-sampling-based
optimisation. It is the de facto default for visualising single-cell
data, though claims that it preserves global structure substantially
better than t-SNE are not well supported; with matched initialisation
(e.g. PCA or Laplacian Eigenmaps), the two methods behave similarly on
global geometry, and both should be interpreted primarily as views of
local structure.

When `use_knn = TRUE` (the default), the kNN graph already stored on the
object (via
[`find_neighbours_sc()`](https://gregorlueg.github.io/bixverse/reference/find_neighbours_sc.md))
is reused, which avoids recomputing nearest neighbours and keeps the
UMAP consistent with any downstream sNN-based clustering. If no kNN is
present, neighbours are computed from the chosen embedding on the fly.

Key parameters to tune: `k` controls the balance between local and
global structure (larger values produce more global layouts), while
`min_dist` and `spread` together control how tightly points are packed
in the embedding. For `MetaCells`, smaller `k` values are often
appropriate given the reduced number of points.

## Usage

``` r
umap_sc(
  object,
  use_knn = TRUE,
  embd_to_use = "pca",
  slot_name = "umap",
  no_embd_to_use = NULL,
  modality = c("rna", "adt", "wnn"),
  n_dim = 2L,
  k = 15L,
  min_dist = 0.5,
  spread = 1,
  knn_method = c("kmknn", "hnsw", "balltree", "annoy", "nndescent", "exhaustive"),
  nn_params = manifoldsR::params_nn(),
  umap_params = manifoldsR::params_umap(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` class.

- use_knn:

  Boolean. Use the kNN graph found in the object. Defaults to `TRUE`. If
  not available, will default to the embedding.

- embd_to_use:

  String. The embedding to use for UMAP. Must be available in the
  object.

- slot_name:

  String. The name of this embedding within the object. Defaults to
  `"umap"`.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

- modality:

  String. On which modality to run the UMAP. One of
  `c("rna", "adt", "wnn")`. The two latter options are only available
  for multi-modal versions with the added data.

- n_dim:

  Integer. Number of UMAP dimensions. Defaults to `2L`.

- k:

  Integer. Number of nearest neighbours. Defaults to `15L`.

- min_dist:

  Numeric. Minimum distance between embedded points. Defaults to `0.5`.

- spread:

  Numeric. Effective scale of embedded points. Defaults to `1.0`.

- knn_method:

  String. Approximate nearest neighbour algorithm. One of `"hnsw"`,
  `"balltree"`, `"annoy"`, `"nndescent"`, or `"exhaustive"`.

- nn_params:

  Named list. See
  [`manifoldsR::params_nn()`](https://gregorlueg.github.io/manifoldsR/reference/params_nn.html).

- umap_params:

  Named list. See
  [`manifoldsR::params_umap()`](https://gregorlueg.github.io/manifoldsR/reference/params_umap.html).

- seed:

  Integer. For reproducibility.

- .verbose:

  Boolean. Controls verbosity.

## Value

The object with a `"umap"` embedding added.
