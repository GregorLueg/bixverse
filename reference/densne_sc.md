# Run den-SNE on a SingleCells/MetaCells object

Wrapper around
[`manifoldsR::densne()`](https://gregorlueg.github.io/manifoldsR/reference/densne.html)
for the `SingleCells` and `MetaCells` classes. den-SNE is t-SNE plus a
density-preserving term: a tight population stays tight in the embedding
and a diffuse one stays diffuse. Plain t-SNE inflates dense clusters and
shrinks sparse ones, so relative cluster sizes on the plot mean nothing;
with den-SNE they do. Setting `lambda = 0` in
[`manifoldsR::params_densne()`](https://gregorlueg.github.io/manifoldsR/reference/params_densne.html)
gives you back plain t-SNE.

Neighbour handling, `approx_type` and the `k` versus `perplexity`
caveats are the same as in
[`tsne_sc()`](https://gregorlueg.github.io/bixverse/reference/tsne_sc.md).

## Usage

``` r
densne_sc(
  object,
  use_knn = FALSE,
  embd_to_use = "pca",
  slot_name = "densne",
  no_embd_to_use = NULL,
  modality = c("rna", "adt", "wnn"),
  n_dim = 2L,
  perplexity = 10,
  approx_type = c("bh", "fft", "fft_3k"),
  knn_method = c("kmknn", "hnsw", "balltree", "annoy", "nndescent", "exhaustive"),
  nn_params = manifoldsR::params_nn(),
  tsne_params = manifoldsR::params_tsne(),
  dens_params = manifoldsR::params_densne(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` class.

- use_knn:

  Boolean. Use the kNN graph found in the object. Defaults to `FALSE`.
  If not available, will default to the embedding.

- embd_to_use:

  String. The embedding to use for den-SNE. Must be available in the
  object.

- slot_name:

  String. The name of this embedding within the object. Defaults to
  `"densne"`.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

- modality:

  String. On which modality to run den-SNE. One of
  `c("rna", "adt", "wnn")`. The two latter options are only available
  for multi-modal versions with the added data.

- n_dim:

  Integer. Number of den-SNE dimensions. Currently only `2L` is
  supported. Defaults to `2L`.

- perplexity:

  Numeric. Perplexity parameter. Typical values between 5 and 50.
  Defaults to `10.0`.

- approx_type:

  String. Approximation method. One of `"bh"` (Barnes-Hut), `"fft"` or
  `"fft_3k"`. Defaults to `"bh"`.

- knn_method:

  String. Approximate nearest neighbour algorithm. One of `"hnsw"`,
  `"balltree"`, `"annoy"`, `"nndescent"`, or `"exhaustive"`.

- nn_params:

  Named list. See
  [`manifoldsR::params_nn()`](https://gregorlueg.github.io/manifoldsR/reference/params_nn.html).

- tsne_params:

  Named list. See
  [`manifoldsR::params_tsne()`](https://gregorlueg.github.io/manifoldsR/reference/params_tsne.html).

- dens_params:

  Named list. The density knobs, see
  [`manifoldsR::params_densne()`](https://gregorlueg.github.io/manifoldsR/reference/params_densne.html).

- seed:

  Integer. For reproducibility.

- .verbose:

  Boolean. Controls verbosity.

## Value

The object with a `"densne"` embedding added.

## References

Narayan, Berger & Cho, Nat. Biotechnol., 2021

## Examples

``` r
# Barnes-Hut den-SNE on the PCA factors
sc <- demo_single_cells()
sc <- densne_sc(sc, .verbose = FALSE)
dim(get_embedding(sc, "densne"))
#> [1] 500   2

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
