# Wrapper function to generate bulk_dge object from h5ad

This is a helper function that can be used to create a `bulk_dge` object
(see [`bulk_dge()`](bulk_dge.md)) directly from h5ad objects.

## Usage

``` r
bulk_dge_from_h5ad(h5_path, .verbose = TRUE)
```

## Arguments

- h5_path:

  String. Path to the h5ad object.

- .verbose:

  Controls verbosity of the function

## Value

`bulk_dge` object.
