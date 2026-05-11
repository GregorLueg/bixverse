# Run t-SNE on a SingleCells/MetaCells object

Wrapper around
[`manifoldsR::tsne()`](https://gregorlueg.github.io/manifoldsR/reference/tsne.html)
for the `SingleCells` and `MetaCells` classes. t-SNE produces a
low-dimensional embedding that emphasises local neighbourhood structure.
Distances between well-separated clusters should not be over-interpreted
quantitatively, but the common claim that t-SNE discards global
structure while UMAP preserves it is largely an artefact of default
initialisations rather than a property of the loss functions themselves.

When `use_knn = FALSE` (the default), the kNN graph already stored on
the object is reused. Otherwise neighbours are computed from the chosen
embedding.

Two approximation strategies are available via `approx_type`: `"bh"`
(Barnes-Hut) is the classical O(n log n) approximation and works well
across a wide range of dataset sizes; `"fft"` (interpolation-based, as
in FIt-SNE) scales better to very large datasets. `perplexity` controls
the bandwidth of the Gaussian kernel used to compute affinities within
the neighbour set (typical values 5-50). When a pre-computed kNN is
supplied via `use_knn = TRUE`, perplexity no longer drives neighbour
retrieval but still shapes the affinity distribution over the retrieved
neighbours; values too close to the kNN size will produce poor results.
With tSNE in particular the rule of thumb is to set k to
`3 * perplexity`. When \`k ≤ perplexity“ the algorithm does not behave
properly anymore, thus, will throw an error.

## Usage

``` r
tsne_sc(
  object,
  use_knn = FALSE,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  n_dim = 2L,
  perplexity = 10,
  approx_type = c("bh", "fft"),
  knn_method = c("kmknn", "hnsw", "balltree", "annoy", "nndescent", "exhaustive"),
  nn_params = manifoldsR::params_nn(),
  tsne_params = manifoldsR::params_tsne(),
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
