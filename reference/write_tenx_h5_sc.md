# Helper function to write data to a 10x CellRanger-style h5 file

Helper function to write data to a 10x CellRanger-style h5 file

## Usage

``` r
write_tenx_h5_sc(
  f_path,
  counts,
  barcodes,
  features,
  version = c("v3", "v2"),
  overwrite = TRUE
)
```

## Arguments

- f_path:

  String. Output path.

- counts:

  Sparse matrix (`dgRMatrix` or `dgCMatrix`), cells x features.

- barcodes:

  Character. Cell barcodes, length `nrow(counts)`.

- features:

  data.table with `id` and `name` of length `ncol(counts)`. For v3 may
  include `feature_type`; defaults to `"Gene Expression"` if absent.

- version:

  One of `"v3"` or `"v2"`.

- overwrite:

  Boolean.

## Value

Invisible.
