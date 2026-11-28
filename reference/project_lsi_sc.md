# Project the query's kept cells onto the reference LSI

Seurat's `ProjectSVD(mode = "lsi")` as
`FindTransferAnchors(reduction = "lsiproject")` runs it: the query's own
TF-IDF times the reference loadings, divided by the singular values,
with no centring or scaling. The reference `"lsi"` embedding is
z-scored, the projection is not; L2 normalise both per cell before
comparing them, as Seurat does. Components the reference dropped are
dropped here too. Either side may be multi-modal or ATAC-only.

## Usage

``` r
project_lsi_sc(reference, query)
```

## Arguments

- reference:

  `SingleCellsAtac` or `SingleCellsMultiModal` class, after
  [`calculate_lsi_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_lsi_sc.md).

- query:

  `SingleCellsAtac` or `SingleCellsMultiModal` class to project, after
  [`run_tfidf_sc()`](https://gregorlueg.github.io/bixverse/reference/run_tfidf_sc.md),
  counted over the reference features.

## Value

The query with the projection in its ATAC cache as `"lsi_projected"`.
