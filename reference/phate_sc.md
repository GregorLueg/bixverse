# Run PHATE on a SingleCells/MetaCells object

Wrapper around
[`manifoldsR::phate()`](https://gregorlueg.github.io/manifoldsR/reference/phate.html)
for the `SingleCells` and `MetaCells` classes. PHATE (Potential of
Heat-diffusion for Affinity-based Trajectory Embedding) produces a
low-dimensional embedding that preserves both local and global structure
by operating on a diffusion process over the data manifold. Unlike UMAP
or t-SNE, PHATE is explicitly designed to reveal continuous progressions
and branching structure, making it the preferred choice for data with
developmental or trajectory-like organisation.

When `use_knn = TRUE` (the default), the kNN graph already stored on the
object is reused; otherwise neighbours are computed from the chosen
embedding. The algorithm then constructs a diffusion operator, raises it
to a power (the diffusion time `t`, see
[`manifoldsR::params_phate()`](https://gregorlueg.github.io/manifoldsR/reference/params_phate.html))
that denoises the manifold, and computes potential distances that are
finally embedded via metric MDS.

Because PHATE inherently smooths over the kNN graph, it pairs naturally
with `MetaCells`: the combination yields a particularly clean view of
continuous biological processes on denoised data.

## Usage

``` r
phate_sc(
  object,
  use_knn = TRUE,
  embd_to_use = "pca",
  slot_name = "phate",
  no_embd_to_use = NULL,
  modality = c("rna", "adt", "wnn"),
  n_dim = 2L,
  k = 5L,
  knn_method = c("kmknn", "hnsw", "balltree", "annoy", "nndescent", "exhaustive"),
  nn_params = manifoldsR::params_nn(),
  phate_params = manifoldsR::params_phate(),
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

  String. The embedding to use for PHATE. Must be available in the
  object.

- slot_name:

  String. The name of this embedding within the object. Defaults to
  `"phate"`.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

- modality:

  String. On which modality to run the UMAP. One of
  `c("rna", "adt", "wnn")`. The two latter options are only available
  for multi-modal versions with the added data.

- n_dim:

  Integer. Number of PHATE dimensions. Currently only `2L` is supported.
  Defaults to `2L`.

- k:

  Integer. Number of nearest neighbours for graph construction. Defaults
  to `5L`.

- knn_method:

  String. Approximate nearest neighbour algorithm. One of `"hnsw"`,
  `"balltree"`, `"annoy"`, `"nndescent"`, or `"exhaustive"`.

- nn_params:

  Named list. See
  [`manifoldsR::params_nn()`](https://gregorlueg.github.io/manifoldsR/reference/params_nn.html).

- phate_params:

  Named list. See
  [`manifoldsR::params_phate()`](https://gregorlueg.github.io/manifoldsR/reference/params_phate.html).

- seed:

  Integer. For reproducibility.

- .verbose:

  Boolean. Controls verbosity.

## Value

The object with a `"phate"` embedding added.
