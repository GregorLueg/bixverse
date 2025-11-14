# Calculate the local pairwise gene-gene correlation

This method implements the HotSpot approach (see DeTomaso, et al.) to
calculate the local gene-gene correlations and their Z-scores.

## Usage

``` r
hotspot_gene_cor_sc(
  object,
  embd_to_use = "pca",
  hotspot_params = params_sc_hotspot(),
  no_embd_to_use = NULL,
  cells_to_take = NULL,
  genes_to_take = NULL,
  streaming = FALSE,
  random_seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- embd_to_use:

  String. The embedding to use. Defaults to `"pca"`.

- hotspot_params:

  List with vision parameters, see
  [`params_sc_hotspot()`](params_sc_hotspot.md) with the following
  elements:

  - model - String. Which of the available models to use for the gene
    expression. Choices are one of `c("danb", "normal", "bernoulli")`.

  - normalise - Boolean. Shall the data be normalised.

  - k - Number of neighbours for the kNN search. Only relevant if you
    set regenerate_knn to `TRUE`.

  - knn_method - String. Which kNN algorithm to use. One of
    `c("annoy", "hnsw", "nndescent")`. Defaults to `"annoy"`. Only
    relevant if you set regenerate_knn to `TRUE`.

  - ann_dist - String. Distance metric for the approximate neighbour
    search. One of `c("cosine", "euclidean")`. Defaults to `"cosine"`.
    Only relevant if you set regenerate_knn to `TRUE`.

  - n_trees - Integer. Number of trees to use for the annoy algorithm.
    Only relevant if you set regenerate_knn to `TRUE`.

  - search_budget - Integer. Search budget per tree for the annoy
    algorithm. Only relevant if you set regenerate_knn to `TRUE`.

  - nn_max_iter - Integer. Maximum iterations for NN Descent. Only
    relevant if you set regenerate_knn to `TRUE` and use `"nndescent"`.

  - rho - Numeric. Sampling rate for NN Descent. Only relevant if you
    set regenerate_knn to `TRUE` and use `"nndescent"`.

  - delta - Numeric. Early termination criterion for NN Descent. Only
    relevant if you set regenerate_knn to `TRUE` and use `"nndescent"`.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

- cells_to_take:

  Optional string vector. If you want to only use selected cells. If
  `NULL` will default to all cells_to_keep in the class.

- genes_to_take:

  Optional string vector. If you wish to limit the search to a subset of
  genes. If `NULL` will default to all genes in the class.

- streaming:

  Boolean. Shall the data be streamed in. Useful for larger data sets.

- random_seed:

  Integer. Used for reproducibility.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A `sc_hotspot` class that can be used for subsequent analysis.
