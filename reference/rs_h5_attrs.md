# Read the attributes of an HDF5 group or dataset into R

Read the attributes of an HDF5 group or dataset into R

## Usage

``` r
rs_h5_attrs(f_path, path)
```

## Arguments

- f_path:

  String. Path to the HDF5 file.

- path:

  String. Path of the group or dataset; `"/"` for the root.

## Value

A named list of vectors. Unsupported dtypes are left out.
