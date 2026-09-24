# Rust implementation of the fgsea::calcGseaStat() function

**\[experimental\]**

## Usage

``` r
rs_calc_gsea_stats(
  stats,
  gs_idx,
  gsea_param,
  return_leading_edge,
  return_all_extremes
)
```

## Arguments

- stats:

  Numeric vector. The gene level statistic. Needs to be sorted in
  descending nature.

- gs_idx:

  Integer vector. The 1-based indices of the gene set genes.

- gsea_param:

  Float. The GSEA parameter. Usually defaults to 1.0.

- return_leading_edge:

  Boolean. Return the leading edge indices.

- return_all_extremes:

  Boolean. Shall the extreme values be returned for plotting.

## Value

List with the following elements

- es Enrichment score for that gene set

- leading_edge 1-based indices of the leading edge genes.

- top Top values of the curve.

- bottom Bottom values of the curve.
