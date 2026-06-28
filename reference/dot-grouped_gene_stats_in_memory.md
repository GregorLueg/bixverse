# Per-group mean expression and percent expressed for in-memory matrices

Per-group mean expression and percent expressed for in-memory matrices

## Usage

``` r
.grouped_gene_stats_in_memory(mat, features, grouping)
```

## Arguments

- mat:

  Numeric matrix. Cells x features.

- features:

  Character vector. Feature ids to extract (pre-matched).

- grouping:

  Factor. Group assignment per cell.

## Value

A long data.table with columns `gene`, `group`, `mean_exp`, `pct_exp`.
