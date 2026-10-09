# Read an h5ad dataset, mapping HDF5 enums to logical

rhdf5 returns HDF5 enums as factors and
[`as.vector()`](https://rdrr.io/r/base/vector.html) on a factor yields
character. anndata only writes enums for booleans.

## Usage

``` r
.read_h5ad_dataset(f_path, ds_path)
```

## Arguments

- f_path:

  String. Path to the h5ad file.

- ds_path:

  String. Full path of the dataset inside the file.

## Value

An atomic vector.
