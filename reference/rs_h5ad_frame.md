# Read an h5ad dataframe group into R

Decodes plain, categorical and nullable columns, old-style
`__categories` and the `_index` attribute or dataset.

## Usage

``` r
rs_h5ad_frame(f_path, group)
```

## Arguments

- f_path:

  String. Path to the h5ad file.

- group:

  String. The dataframe group, e.g. `"obs"`.

## Value

A list with `cols` (named list of vectors, categoricals as factors, in
HDF5 name order with column groups first), `idx` (the index, or `NULL`)
and `skipped` (columns with an unsupported encoding).
