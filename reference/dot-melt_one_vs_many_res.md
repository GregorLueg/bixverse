# Melt one reference arm of the one-vs-many marker results

Takes the raw Rust output of
[`rs_calculate_dge_one_vs_many()`](https://gregorlueg.github.io/bixverse/reference/rs_calculate_dge_one_vs_many.md)
for a single reference group and splits it into the per-rival statistics
and the per-gene summaries across those rivals. The per-comparison
elements come back comparison-major, so the gene names are recycled and
each rival name covers one block of genes.

## Usage

``` r
.melt_one_vs_many_res(rs_res, gene_names, ref_grp, rival_grps)
```

## Arguments

- rs_res:

  List. The raw return of
  [`rs_calculate_dge_one_vs_many()`](https://gregorlueg.github.io/bixverse/reference/rs_calculate_dge_one_vs_many.md).
  Must have at least one gene that passed the proportion filter.

- gene_names:

  Character vector. All gene names in the original gene order, subset
  with `rs_res$genes_to_keep`.

- ref_grp:

  String. Name of the reference group.

- rival_grps:

  Character vector. Names of the comparison groups, in the order they
  were handed to Rust.

## Value

A list with elements `summary` and `per_comparison`, both data.tables.
