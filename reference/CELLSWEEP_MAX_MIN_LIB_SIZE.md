# Largest smallest-library-size that still looks like raw ingest

A raw, unfiltered barcode list always contains near-empty droplets. If
the minimum library size across the whole object is above this, a QC
cutoff was applied at load time and the empty droplets are gone.

## Usage

``` r
CELLSWEEP_MAX_MIN_LIB_SIZE
```
