# Column the HVG statistics are ranked on

Maps an HVG method to the column of the `rs_sc_hvg` / `rs_mc_hvg` output
that genes are ranked on, highest first.

## Usage

``` r
hvg_rank_col(hvg_method)
```

## Arguments

- hvg_method:

  String. One of `c("vst", "dispersion", "meanvarbin", "scran")`.

## Value

The name of the ranking column.
