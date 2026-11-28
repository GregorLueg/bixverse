# The store a marker test reads

The store a marker test reads

## Usage

``` r
.de_store(object, modality)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- modality:

  String. One of `c("rna", "atac")`.

## Value

A list with `cells` and `genes` (the store's two files), `names` (the
feature names in store order), `id_col` (`"gene_id"` or `"feature_id"`)
and `atac` (whether the store holds ATAC features).
