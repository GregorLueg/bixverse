# Total counts per feature and the top features

**\[experimental\]** Signac's `FindTopFeatures()`.

## Usage

``` r
rs_atac_feature_totals(f_path_genes, cells, min_cutoff)
```

## Arguments

- f_path_genes:

  String. Gene-based file.

- cells:

  Optional integer vector. Cells to sum over (0-based), `NULL` for all.

- min_cutoff:

  Optional string. `"q<N>"` keeps features above the N-th percentile, a
  number those with more counts, `NULL` or `"all"` every feature.

## Value

A list with `count`, `percentile` and `top_feature`.
