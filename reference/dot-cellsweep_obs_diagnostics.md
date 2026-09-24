# Assemble the per-barcode CellSweep diagnostics

Stitches the per-sample fits back into one table in output order.
Per-sample scalars are broadcast across that sample's barcodes, which
keeps them queryable in obs without needing a second table.

## Usage

``` r
.cellsweep_obs_diagnostics(res, celltype_levels)
```

## Arguments

- res:

  List. The Rust return value.

- celltype_levels:

  Character vector. Factor levels the `z_hat` codes index into.

## Value

A data.table with one row per written barcode.
