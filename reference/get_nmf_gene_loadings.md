# Get the NMF gene loadings

Getter function to extract the gene loadings matrix (features x k) from
a bulk NMF fit stored in
[`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).
If no NMF fit is present, returns `NULL` with a warning.

## Usage

``` r
get_nmf_gene_loadings(object)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

## Value

A features x k numeric matrix (if found) or `NULL`.
