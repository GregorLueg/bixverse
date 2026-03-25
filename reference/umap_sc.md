# Run UMAP on a SingleCells object

Wrapper around
[`manifoldsR::umap()`](https://gregorlueg.github.io/manifoldsR/reference/umap.html)
for the `SingleCells` class.

## Usage

``` r
umap_sc(
  object,
  use_knn = TRUE,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  n_dim = 2L,
  k = 15L,
  min_dist = 0.5,
  spread = 1,
  knn_method = c("hnsw", "balltree", "annoy", "nndescent", "exhaustive"),
  nn_params = manifoldsR::params_nn(),
  umap_params = manifoldsR::params_umap(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- use_knn:

  Boolean. Use the kNN graph found in the object. Defaults to `TRUE`. If
  not available, will default to the embedding.

- embd_to_use:

  String. The embedding to use for UMAP. Must be available in the
  object.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

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
