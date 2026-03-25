# Resolve the index column from an h5ad group

Checks for the `_index` attribute first (new anndata spec), then falls
back to the `_index` dataset (old spec). Returns NULL if neither exists.

## Usage

``` r
.resolve_h5_index(f_path, group_path, h5_content)
```

## Arguments

- f_path:

  Path to the h5ad file.

- group_path:

  The HDF5 group (e.g. "/obs", "/var").

- h5_content:

  data.table from rhdf5::h5ls().

## Value

A list with `idx` (character vector or NULL) and `idx_col` (the dataset
name used, or NULL).
