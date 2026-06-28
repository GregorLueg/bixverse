# Calculate the local pairwise gene-gene correlation

This method implements the HotSpot approach (see DeTomaso, et al.) to
calculate the local gene-gene correlations and their Z-scores.

## Usage

``` r
hotspot_gene_cor_sc(
  object,
  embd_to_use = "pca",
  use_knn = TRUE,
  hotspot_params = params_sc_hotspot(),
  no_embd_to_use = NULL,
  cells_to_take = NULL,
  genes_to_take = NULL,
  streaming = NULL,
  working_mem_gb = 4,
  random_seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- embd_to_use:

  String. The embedding to use. Defaults to `"pca"`.

- use_knn:

  Boolean. Shall the internal kNN be used. If set to yes, you need to
  ensure consistency. If you provide `cells_to_take`, the function will
  regenerate the kNN graph with these cells.

- hotspot_params:

  List with hotspot parameters, see
  [`params_sc_hotspot()`](https://gregorlueg.github.io/bixverse/reference/params_sc_hotspot.md)
  with the following elements:

  - model - String. Which of the available models to use for the gene
    expression. Choices are one of `c("danb", "normal", "bernoulli")`.

  - normalise - Boolean. Shall the data be normalised.

  - knn - List of kNN parameters. See
    [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
    for available parameters and their defaults.

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

  Optional Boolean. Shall the data be streamed in. Useful for larger
  data sets where you wish to avoid loading in the whole data. If
  `NULL`, will automatically detect.

- working_mem_gb:

  Numeric. Approximate working memory (GB) the streaming pair path may
  use for resident gene panels. Ignored when `streaming` is `FALSE`.
  Larger values mean fewer disk re-reads. Note this excludes the two
  dense N_genes x N_genes output matrices, which scale with
  `genes_to_use`. Defaults to `4` (4 GB of memory allocated).

- random_seed:

  Integer. Used for reproducibility.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

A `sc_hotspot` class that can be used for subsequent analysis.

## Details

Should a gene not be found in sufficient cells, the pairs with this gene
will be set to 0. Please ensure prior to running the function that you
are only calculating gene-gene auto-correlations that occur in
sufficient cells.
