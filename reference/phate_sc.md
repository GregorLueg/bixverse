# Run PHATE on a SingleCells object

Wrapper around
[`manifoldsR::phate()`](https://gregorlueg.github.io/manifoldsR/reference/phate.html)
for the `SingleCells` class.

## Usage

``` r
phate_sc(
  object,
  use_knn = TRUE,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  n_dim = 2L,
  k = 5L,
  knn_method = c("hnsw", "balltree", "annoy", "nndescent", "exhaustive"),
  nn_params = manifoldsR::params_nn(),
  phate_params = manifoldsR::params_phate(),
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

  String. The embedding to use for PHATE. Must be available in the
  object.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

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
