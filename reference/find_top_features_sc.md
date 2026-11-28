# Find the most accessible features

Signac's `FindTopFeatures()`. Totals are summed over the kept cells; the
ATAC feature table gains `count`, `percentile` and `top_feature`, which
[`calculate_lsi_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_lsi_sc.md)
then uses.

## Usage

``` r
find_top_features_sc(object, min_cutoff = "q5")
```

## Arguments

- object:

  `SingleCellsAtac` or `SingleCellsMultiModal` class with ATAC.

- min_cutoff:

  `"qN"` keeps features above the N-th percentile of the totals, a
  number those with more counts, `NULL` or `"all"` all.

## Value

The object with the columns added to the ATAC feature table.
