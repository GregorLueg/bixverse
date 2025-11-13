# Helper function to get the dimensions and compressed sparse format

Helper function to get the dimensions and compressed sparse format

## Usage

``` r
get_h5ad_dimensions(f_path)
```

## Arguments

- f_path:

  File path to the `.h5ad` file.

## Value

A list with the following elements:

- dims - Dimensions of the stored data in the h5ad file.

- type - Was the data stored in CSR (indptr = cells) or CSC (indptr =
  genes).
