# Read one HDF5 dataset into R

1-D datasets and encoded columns (categorical, nullable) come back as
vectors. A 2-D numeric dataset comes back as a double matrix with the
dimensions reversed (the on-disk row-major `rows x cols` read as R's
column-major `cols x rows`).

## Usage

``` r
rs_h5_read(f_path, path)
```

## Arguments

- f_path:

  String. Path to the HDF5 file.

- path:

  String. Path of the dataset inside the file.

## Value

The vector or matrix.
