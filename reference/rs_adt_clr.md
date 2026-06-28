# Applies CLR normalisation on ADT counts

**\[experimental\]**

## Usage

``` r
rs_adt_clr(counts, seurat_clr)
```

## Arguments

- counts:

  R matrix of shape cells x features.

- seurat_clr:

  Logical; if TRUE uses the Seurat variant (non-negative), if FALSE uses
  proper CLR (mean-centred log, can be negative).

## Value

CLR-transformed matrix.
