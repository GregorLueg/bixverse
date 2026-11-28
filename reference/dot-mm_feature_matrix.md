# Kept-cell matrix of ADT or ATAC features of a multi-modal object

The normalised ADT counts, or for ATAC the normalised counts of the
requested features read from the ATAC store (log-CPM, or TF-IDF after
[`run_tfidf_sc()`](https://gregorlueg.github.io/bixverse/reference/run_tfidf_sc.md)),
as a dense matrix over the kept cells.

## Usage

``` r
.mm_feature_matrix(object, modality, features)
```

## Arguments

- object:

  `SingleCellsMultiModal` class.

- modality:

  String. One of `c("adt", "atac")`.

- features:

  Character vector. The requested feature ids.

## Value

A list with `mat` (kept cells x features, cell ids as row names, feature
ids as column names) and `features` (the matched feature ids).
