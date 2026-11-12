# List the members of an HDF5 group

List the members of an HDF5 group

## Usage

``` r
rs_h5_ls(f_path, group)
```

## Arguments

- f_path:

  String. Path to the HDF5 file.

- group:

  String. Path of the group; `"/"` for the root.

## Value

A list with `name`, `otype` (`"group"` or `"dataset"`) and `rank`
(number of dimensions, 0 for groups and scalars).
