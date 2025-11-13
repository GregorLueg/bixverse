# Wrapper function to provide data for mtx-based loading

Wrapper function to provide data for mtx-based loading

## Usage

``` r
params_sc_mtx_io(path_mtx, path_obs, path_var, cells_as_rows, has_hdr)
```

## Arguments

- path_mtx:

  String. Path to the .mtx file

- path_obs:

  String. Path to the file containing cell/barcode info.

- path_var:

  String. Path to the file containing gene/variable info.

- cells_as_rows:

  Boolean. Do cells represent the rows or columns.

- has_hdr:

  Boolean. Do the plain text files have headers.

## Value

A list with the mtx loading parameters for usage in subsequent
functions.
