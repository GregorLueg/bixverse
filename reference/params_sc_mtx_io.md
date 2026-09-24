# Wrapper function to provide data for mtx-based loading

Wrapper function to provide data for mtx-based loading

## Usage

``` r
params_sc_mtx_io(path_mtx, path_obs, path_var, cells_as_rows, has_hdr)
```

## Arguments

- path_mtx:

  Any. Path to the .mtx file Required.

- path_obs:

  Any. Path to the file containing cell/barcode info. Required.

- path_var:

  Any. Path to the file containing gene/variable info. Required.

- cells_as_rows:

  Boolean. Do cells represent the rows or columns. Required.

- has_hdr:

  Boolean. Do the plain text files have headers. Required.

## Value

A named list with the following elements:

- path_mtx - Any. Path to the .mtx file Required.

- path_obs - Any. Path to the file containing cell/barcode info.
  Required.

- path_var - Any. Path to the file containing gene/variable info.
  Required.

- cells_as_rows - Boolean. Do cells represent the rows or columns.
  Required.

- has_hdr - Boolean. Do the plain text files have headers. Required.
