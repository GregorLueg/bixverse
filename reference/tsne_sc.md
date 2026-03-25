# Run t-SNE on a SingleCells object

Wrapper around
[`manifoldsR::tsne()`](https://gregorlueg.github.io/manifoldsR/reference/tsne.html)
for the `SingleCells` class.

## Usage

``` r
tsne_sc(
  object,
  use_knn = TRUE,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  n_dim = 2L,
  perplexity = 30,
  approx_type = c("bh", "fft"),
  knn_method = c("hnsw", "balltree", "annoy", "nndescent", "exhaustive"),
  nn_params = manifoldsR::params_nn(),
  tsne_params = manifoldsR::params_tsne(),
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

  String. The embedding to use for t-SNE. Must be available in the
  object.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

- n_dim:

  Integer. Number of t-SNE dimensions. Currently only `2L` is supported.
  Defaults to `2L`.

- perplexity:

  Numeric. Perplexity parameter. Typical values between 5 and 50.
  Defaults to `30.0`.

- approx_type:

  String. Approximation method. One of `"bh"` (Barnes-Hut) or `"fft"`.
  Defaults to `"bh"`.

- knn_method:

  String. Approximate nearest neighbour algorithm. One of `"hnsw"`,
  `"balltree"`, `"annoy"`, `"nndescent"`, or `"exhaustive"`.

- nn_params:

  Named list. See
  [`manifoldsR::params_nn()`](https://gregorlueg.github.io/manifoldsR/reference/params_nn.html).

- tsne_params:

  Named list. See
  [`manifoldsR::params_tsne()`](https://gregorlueg.github.io/manifoldsR/reference/params_tsne.html).

- seed:

  Integer. For reproducibility.

- .verbose:

  Boolean. Controls verbosity.

## Value

The object with a `"tsne"` embedding added.
