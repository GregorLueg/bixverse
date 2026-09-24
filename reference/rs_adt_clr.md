# Applies CLR normalisation on ADT counts

**\[experimental\]**

Normalises each cell (row) separately.

## Usage

``` r
rs_adt_clr(counts, seurat_clr)
```

## Arguments

- counts:

  Numerical matrix of shape cells x features.

- seurat_clr:

  Boolean. If `TRUE` uses the Seurat variant `log1p(x / g)`
  (non-negative); if `FALSE` uses the proper CLR
  `log1p(x) - mean(log1p(x))` (mean-centred, can be negative).

## Value

Numerical matrix of cells x features with the CLR-transformed values.
