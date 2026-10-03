# Melt the one-vs-many marker results

Takes the raw Rust output of
[`rs_calculate_dge_one_vs_many()`](https://gregorlueg.github.io/bixverse/reference/rs_calculate_dge_one_vs_many.md)
and splits it into the per-rival statistics and the per-gene summaries
across those rivals. Both come back with 0-indexed group indices and are
flattened reference-major, so the kept gene names are recycled per
block.

## Usage

``` r
.melt_one_vs_many_res(rs_res, gene_names, grp_names)
```

## Arguments

- rs_res:

  List. The raw return of
  [`rs_calculate_dge_one_vs_many()`](https://gregorlueg.github.io/bixverse/reference/rs_calculate_dge_one_vs_many.md).
  Must have at least one gene that passed the proportion filter.

- gene_names:

  Character vector. All gene names in the original gene order, subset
  with `rs_res$genes_to_keep`.

- grp_names:

  Character vector. Names of all groups, in the order they were handed
  to Rust.

## Value

A list with elements `summary` and `per_comparison`, both data.tables.
