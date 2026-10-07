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
    `c("vst", "meanvarbin", "dispersion", "scran", "residual")`

  - loess_span - The span for the loess function to standardise the
    variance (`"vst"`), or of the lowess trend (`"scran"`)

  - num_bin - Integer. Not yet implemented.

  - bin_method - String. One of `c("equal_width", "equal_freq")`. Not
    implemented yet.

  - mean_filter, min_mean, transform, use_min_width, min_width,
    min_window_count - The `"scran"` trend parameters, see
    [`params_hvg_scran_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_hvg_scran_defaults.md)

- streaming:

  Optional Boolean. Shall the data be streamed in. Useful for larger
  data sets where you wish to avoid loading in the whole data. If
  `NULL`, will automatically detect. Not used for `MetaCells`.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

It will add the per-gene HVG statistics to the var table: `mean`, `var`,
`var_exp` and `var_std` for `"vst"`; `mean`, `dispersion`,
`dispersion_scaled` and `bin` for `"meanvarbin"` and `"dispersion"`;
`scran_mean`, `scran_var`, `scran_fitted` and `scran_residual` (log2
scale) for `"scran"`.

## Examples

``` r
# the twenty most variable genes by the vst method
sc <- demo_single_cells(prepped = FALSE)
sc <- find_hvg_sc(sc, hvg_no = 20L, .verbose = FALSE)
length(get_hvg(sc))
#> [1] 20

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
