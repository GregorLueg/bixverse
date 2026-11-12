# Version and dimensions of a 10x CellRanger h5 file

Version and dimensions of a 10x CellRanger h5 file

## Usage

``` r
rs_tenx_h5_metadata(f_path)
```

## Arguments

- f_path:

  String. Path to the 10x h5 file.

## Value

A list with `version` (`"v2"` or `"v3"`), `n_genes` and `n_cells`.
