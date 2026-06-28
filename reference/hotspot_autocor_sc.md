# Calculate the local auto-correlation of a gene

This method implements the HotSpot approach (see DeTomaso, et al.) to
calculate the auto-correlation of a given gene in the kNN graph based on
the chosen embedding. This can be used to identify genes that have
strong local correlations and vary across the kNN graph.

## Usage

``` r
hotspot_autocor_sc(
  object,
  embd_to_use = "pca",
  use_knn = TRUE,
  hotspot_params = params_sc_hotspot(),
  no_embd_to_use = NULL,
  cells_to_take = NULL,
  genes_to_take = NULL,
  streaming = NULL,
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

- random_seed:

  Integer. Used for reproducibility.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

A data.table with the auto-correlations on a per gene basis and various
statistics.

## Details

Should a gene not be found in sufficient cells, the gene will be
automatically filtered out from the results. This can occur for example
if you have filtered out the cells that contain a given gene. The
underlying genes are still available, but the cells that might contain
them are not included.
