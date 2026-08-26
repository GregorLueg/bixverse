# Coerce a document-term matrix into the sparse list Rust wants

Accepts a `dgCMatrix`, a `dgRMatrix` or a dense numeric/logical matrix.
Logical input is what
[`binarise_regulon_activity()`](https://gregorlueg.github.io/bixverse/reference/binarise_regulon_activity.md)
hands back, so it is promoted rather than rejected.

## Usage

``` r
.lda_matrix_to_list(x)
```

## Arguments

- x:

  The document-term matrix.

## Value

A list with `data`, `indices`, `indptr`, `cs_type`, `nrow` and `ncol`,
plus the row and column names as attributes `doc_ids` and `term_ids`.
