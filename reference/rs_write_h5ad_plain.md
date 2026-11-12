# Write a minimal h5ad file

Writes `X` plus `obs` / `var` as plain datasets without anndata encoding
attributes. The first column of `obs` and `var` becomes `_index`.

## Usage

``` r
rs_write_h5ad_plain(
  f_path,
  x_type,
  data,
  indices,
  indptr,
  n_obs,
  n_var,
  obs,
  var
)
```

## Arguments

- f_path:

  String. Output file, overwritten.

- x_type:

  String. One of `c("csr", "csc", "dense")`.

- data:

  Numeric. Non-zero values, or the dense values cells x genes in
  row-major order.

- indices:

  Integer. 0-based minor indices; ignored for dense.

- indptr:

  Integer. 0-based major offsets; ignored for dense.

- n_obs:

  Integer. Number of cells.

- n_var:

  Integer. Number of genes.

- obs:

  Named list of cell columns, identifier first.

- var:

  Named list of gene columns, identifier first.
