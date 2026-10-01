# Run densMAP on a SingleCells/MetaCells object

Wrapper around
[`manifoldsR::densmap()`](https://gregorlueg.github.io/manifoldsR/reference/densmap.html)
for the `SingleCells` and `MetaCells` classes. densMAP is UMAP plus a
density-preserving term: a tight population stays tight in the embedding
and a diffuse one stays diffuse. With plain UMAP the relative size of a
cluster on the plot tells you nothing, with densMAP it does. Setting
`lambda = 0` in
[`manifoldsR::params_densmap()`](https://gregorlueg.github.io/manifoldsR/reference/params_densmap.html)
gives you back plain UMAP.

Neighbour handling is the same as in
[`umap_sc()`](https://gregorlueg.github.io/bixverse/reference/umap_sc.md):
with `use_knn = TRUE` (the default) the cached kNN graph is reused,
otherwise neighbours are computed from the chosen embedding.

## Usage

``` r
densmap_sc(
  object,
  use_knn = TRUE,
  embd_to_use = "pca",
  slot_name = "densmap",
  no_embd_to_use = NULL,
  modality = c("rna", "adt", "wnn"),
  n_dim = 2L,
  k = 15L,
  min_dist = 0.5,
  spread = 1,
  knn_method = c("kmknn", "hnsw", "annoy", "nndescent", "balltree", "ivf", "exhaustive"),
  nn_params = manifoldsR::params_nn(),
  umap_params = manifoldsR::params_umap(),
  dens_params = manifoldsR::params_densmap(),
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

  String. The embedding to use for densMAP. Must be available in the
  object.

- slot_name:

  String. The name of this embedding within the object. Defaults to
  `"densmap"`.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

- modality:

  String. On which modality to run densMAP. One of
  `c("rna", "adt", "wnn")`. The two latter options are only available
  for multi-modal versions with the added data.

- n_dim:

  Integer. Number of densMAP dimensions. Defaults to `2L`.

- k:

  Integer. Number of nearest neighbours. Defaults to `15L`.

- min_dist:

  Numeric. Minimum distance between embedded points. Defaults to `0.5`.

- spread:

  Numeric. Effective scale of embedded points. Defaults to `1.0`.

- knn_method:

  String. Approximate nearest neighbour algorithm. One of
  `c("kmknn", "hnsw", "annoy", "nndescent", "balltree", "ivf", "exhaustive")`.
  Defaults to `"kmknn"`. Only used when neighbours are computed from the
  embedding.

- nn_params:

  Named list. See
  [`manifoldsR::params_nn()`](https://gregorlueg.github.io/manifoldsR/reference/params_nn.html).

- umap_params:

  Named list. See
  [`manifoldsR::params_umap()`](https://gregorlueg.github.io/manifoldsR/reference/params_umap.html).

- dens_params:

  Named list. The density knobs, see
  [`manifoldsR::params_densmap()`](https://gregorlueg.github.io/manifoldsR/reference/params_densmap.html).

- seed:

  Integer. For reproducibility.

- .verbose:

  Boolean. Controls verbosity.

## Value

The object with a `"densmap"` embedding added.

## References

Narayan, Berger & Cho, Nat. Biotechnol., 2021

## Examples

``` r
# densMAP off the cached kNN graph
sc <- demo_single_cells()
sc <- densmap_sc(sc, .verbose = FALSE)
dim(get_embedding(sc, "densmap"))
#> [1] 500   2

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
