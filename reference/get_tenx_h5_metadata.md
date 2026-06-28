# Detect 10x h5 version and read dimensions

Detect 10x h5 version and read dimensions

## Usage

``` r
get_tenx_h5_metadata(f_path)
```

## Arguments

- f_path:

  Path to the 10x CellRanger h5 file.

## Value

A list with `version` (`"v2"`/`"v3"`), `n_cells` and `n_genes`.
