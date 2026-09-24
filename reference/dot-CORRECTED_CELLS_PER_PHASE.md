# Cells held in memory per phase when transposing a corrected store

Peak memory is roughly `cells_per_phase * mean_genes_per_cell * 12`
bytes, so 50k cells is a few hundred megabytes at a typical density. The
source is re-read once per phase, which is the trade.

## Usage

``` r
.CORRECTED_CELLS_PER_PHASE
```
