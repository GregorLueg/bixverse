# Read obs and var tables and metadata from an h5ad file

Useful for exploring the data stored in an h5ad file.

## Usage

``` r
read_h5ad_metadata(f_path)
```

## Arguments

- f_path:

  File path to the `.h5ad` file.

## Value

A list with:

- obs - data.table of cell-level metadata

- var - data.table of gene-level metadata

- dims - named integer vector c(obs, var)

- type - "CSR" or "CSC"
