# Identify HVGs

This is a helper function to identify highly variable genes. At the
moment the implementation has only the VST-based version (known as
Seurat v3). The other methods will be implemented in the future.

## Usage

``` r
find_hvg_sc(
  object,
  hvg_no = 2000L,
  hvg_params = params_sc_hvg(),
  streaming = FALSE,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- hvg_no:

  Integer. Number of highly variable genes to include. Defaults to
  `2000L`.

- hvg_params:

  List, see [`params_sc_hvg()`](params_sc_hvg.md). This list contains

  - method - Which method to use. One of
    `c("vst", "meanvarbin", "dispersion")`

  - loess_span - The span for the loess function to standardise the
    variance

  - num_bin - Integer. Not yet implemented.

  - bin_method - String. One of `c("equal_width", "equal_freq")`. Not
    implemented yet.

- streaming:

  Boolean. Shall the genes be streamed in. Useful for larger data sets
  where you wish to avoid loading in the whole data. Defaults to
  `FALSE`.

- .verbose:

  Boolean. Controls verbosity and returns run times.

## Value

It will add the mean, var, var_exp, var_std of each gene to the the var
table.
