# Read the columns of an h5ad dataframe group

Columns stored as groups come first, then plain datasets, as h5ls lists
them. Old-style categoricals (codes in a dataset, levels under
`__categories`) are resolved. The index column is returned separately.

## Usage

``` r
.read_h5ad_frame(f_path, group_path, h5_content)
```

## Arguments

- f_path:

  String. Path to the h5ad file.

- group_path:

  String. The dataframe group, e.g. `"/obs"`.

- h5_content:

  data.table. Output of
  [`rhdf5::h5ls()`](https://huber-group-embl.github.io/rhdf5/reference/h5ls.html)
  on `f_path`.

## Value

A list with `cols` (named list of column vectors) and `idx` (character
vector or `NULL`).
