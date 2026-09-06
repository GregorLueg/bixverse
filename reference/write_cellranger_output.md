# Helper function to write data to a cell ranger like output

This is a helper to write synthetic data to cell ranger like output,
i.e., an .mtx file, an barcodes.csv (or .tsv) and a features.csv (or
.tsv).

## Usage

``` r
write_cellranger_output(
  f_path,
  counts,
  obs,
  var,
  format_type = c("csv", "tsv"),
  rows = c("cells", "genes"),
  overwrite = TRUE,
  .verbose = TRUE
)
```

## Arguments

- f_path:

  String. The filepath to which to save the data

- counts:

  Sparse matrix. Needs to be of class `dgRMatrix` or `dgCMatrix`.

- obs:

  data.table. The observations. Needs to have
  `nrow(obs) == nrow(counts)`.

- var:

  data.table. The variable data. Needs to have
  `ncol(var) == ncol(counts)`.

- format_type:

  String. One of `c("csv", "tsv")`. Shall the data be saved in TSV or
  CSV.

- rows:

  String. One of `c("cells", "genes")`. Shall the rows represent cells
  or genes in the .mtx file.

- overwrite:

  Boolean. Shall any found h5ad file be overwritten.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

Returns invisible

## Examples

``` r
# the 10x trio: an .mtx plus barcode and feature tables
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
)
dir_out <- tempfile("cellranger")
dir.create(dir_out, recursive = TRUE)
write_cellranger_output(
  f_path = dir_out,
  counts = data$counts,
  obs = data$obs,
  var = data$var,
  rows = "cells",
  format_type = "csv",
  .verbose = FALSE
)
list.files(dir_out)
#> [1] "barcodes.csv" "features.csv" "matrix.mtx"  

unlink(dir_out, recursive = TRUE, force = TRUE)
```
