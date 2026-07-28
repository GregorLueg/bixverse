# Get the stabilised NMF diagnostics

Getter function to extract per-run losses, convergence flags, and the
best-run index from a
[`stabilised_nmf_bulk()`](https://gregorlueg.github.io/bixverse/reference/stabilised_nmf_bulk.md)
fit. Returns `NULL` for a single-run fit.

## Usage

``` r
get_nmf_stability(object)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

## Value

A list with `losses`, `converged`, `best_idx` (if found) or `NULL`.
