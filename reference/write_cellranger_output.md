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
