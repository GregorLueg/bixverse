# Scatter per-gene values into a full-length vector

The variable table takes one row per gene in the store, but a residual
model only covers the genes it retained. Everything else is `NA` rather
than zero, which would read as a real measurement of no variance.

## Usage

``` r
.scatter_gene_values(values, gene_indices, n_genes)
```

## Arguments

- values:

  Numeric. One value per entry in `gene_indices`.

- gene_indices:

  Integer. The 0-based genes `values` belongs to.

- n_genes:

  Integer. Total number of genes in the store.

## Value

A numeric vector of length `n_genes`.
