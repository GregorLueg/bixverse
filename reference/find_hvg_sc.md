# Identify HVGs

This is a helper function to identify highly variable genes for
`SingleCells` (using the Rust-based streaming of data) or `MetaCells`.

## Usage

``` r
find_hvg_sc(
  object,
  hvg_no = 2000L,
  hvg_params = params_sc_hvg(),
  streaming = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` (or potentially other) class.

- hvg_no:

  Integer. Number of highly variable genes to include. Defaults to
  `2000L`.

- hvg_params:

  List, see
  [`params_sc_hvg()`](https://gregorlueg.github.io/bixverse/reference/params_sc_hvg.md).
  This list contains

  - method - Which method to use. One of
    `c("vst", "meanvarbin", "dispersion")`

  - loess_span - The span for the loess function to standardise the
    variance

  - num_bin - Integer. Not yet implemented.

  - bin_method - String. One of `c("equal_width", "equal_freq")`. Not
    implemented yet.

- streaming:

  Optional Boolean. Shall the data be streamed in. Useful for larger
  data sets where you wish to avoid loading in the whole data. If
  `NULL`, will automatically detect. Not used for `MetaCells`.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

It will add the mean, var, var_exp, var_std of each gene to the the var
table.
