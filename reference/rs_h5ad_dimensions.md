# Dimensions and storage layout of an h5ad file

**\[experimental\]** Reads the obs and var lengths and how `X` is
stored, without touching the counts. Sparse `X` takes its layout from
the `encoding-type` attribute, dense `X` from comparing its shape
against the obs and var lengths.

## Usage

``` r
rs_h5ad_dimensions(f_path)
```

## Arguments

- f_path:

  String. Path to the h5ad file.

## Value

A list with

- n_obs - Integer. Number of cells (rows of obs).

- n_var - Integer. Number of genes (rows of var).

- format - String. One of `c("CSR", "CSC", "DENSE_ROW", "DENSE_COL")`.
