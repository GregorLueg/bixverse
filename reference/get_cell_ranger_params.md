# Helper to generate cell ranger input parameters

Resolves the three files a Cell Ranger MTX directory holds and wraps
them into the parameter list
[`load_mtx()`](https://gregorlueg.github.io/bixverse/reference/load_mtx.md)
wants. Handles the v2 naming (`genes.tsv`) and the v3 one
(`features.tsv`), and the `.csv` variant that
[`write_cellranger_output()`](https://gregorlueg.github.io/bixverse/reference/write_cellranger_output.md)
can emit.

## Usage

``` r
get_cell_ranger_params(dir_data, cells_as_rows = FALSE, has_hdr = FALSE)
```

## Arguments

- dir_data:

  String. The directory with the Cell Ranger outputs

- cells_as_rows:

  Boolean. Are the cells the rows of the matrix? Cell Ranger writes
  genes x cells, so this defaults to `FALSE`. Set to `TRUE` for output
  of
  [`write_cellranger_output()`](https://gregorlueg.github.io/bixverse/reference/write_cellranger_output.md)
  written with `rows = "cells"`.

- has_hdr:

  Boolean. Do the barcode and feature files carry a header row? Cell
  Ranger writes none, so this defaults to `FALSE`.
  [`write_cellranger_output()`](https://gregorlueg.github.io/bixverse/reference/write_cellranger_output.md)
  does write one.

## Value

A list based on
[`params_sc_mtx_io()`](https://gregorlueg.github.io/bixverse/reference/params_sc_mtx_io.md).

## Examples

``` r
# round trip through the package's own writer
dir <- tempfile("cellranger")
dir.create(dir)
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 50L, n_genes = 40L)
)
write_cellranger_output(
  f_path = dir,
  counts = data$counts,
  obs = data$obs,
  var = data$var,
  rows = "genes",
  format_type = "tsv",
  .verbose = FALSE
)
str(get_cell_ranger_params(dir, has_hdr = TRUE))
#> List of 5
#>  $ path_mtx     : chr "/tmp/RtmpvYy6aB/cellranger47b12f5fe1ba/matrix.mtx"
#>  $ path_obs     : chr "/tmp/RtmpvYy6aB/cellranger47b12f5fe1ba/barcodes.tsv"
#>  $ path_var     : chr "/tmp/RtmpvYy6aB/cellranger47b12f5fe1ba/features.tsv"
#>  $ cells_as_rows: logi FALSE
#>  $ has_hdr      : logi TRUE

unlink(dir, recursive = TRUE, force = TRUE)
```
