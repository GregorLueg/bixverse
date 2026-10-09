# Read one column of an h5ad dataframe

Handles plain arrays, categoricals (`categories` + `codes`) and the
nullable encodings (`values` + `mask`), which anndata \>= 0.13 also uses
for every non-categorical string column.

## Usage

``` r
.read_h5ad_column(f_path, col_path, h5_content)
```

## Arguments

- f_path:

  String. Path to the h5ad file.

- col_path:

  String. Full path of the column inside the file.

- h5_content:

  data.table. Output of
  [`rhdf5::h5ls()`](https://huber-group-embl.github.io/rhdf5/reference/h5ls.html)
  on `f_path`.

## Value

An atomic vector or factor, or `NULL` for an unsupported encoding.
