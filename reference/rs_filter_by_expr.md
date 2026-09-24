# Filter lowly expressed genes

**\[experimental\]** edgeR's `filterByExpr`, via the `edge-rs` crate.

## Usage

``` r
rs_filter_by_expr(
  counts,
  group,
  lib_size,
  min_count,
  min_total_count,
  min_prop
)
```

## Arguments

- counts:

  Integer or double matrix. Raw counts of genes x samples.

- group:

  Integer vector or NULL. 1-based group code per sample (e.g.
  `as.integer(factor(x))`). If NULL, all samples form one group.

- lib_size:

  Numeric vector or NULL. Library size per sample. NULL uses the column
  sums.

- min_count:

  Numeric. Minimum count in the median-sized library.

- min_total_count:

  Numeric. Minimum total count across all samples.

- min_prop:

  Numeric. Fraction of the smallest group size beyond edgeR's `large.n`
  that still has to express the gene.

## Value

Boolean vector, one per gene. `TRUE` if the gene is kept.

## References

Chen, Lun and Smyth, F1000Research, 2016
