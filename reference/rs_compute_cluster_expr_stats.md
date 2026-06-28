# Compute cluster statistics for NicheNet prioritisation

**\[experimental\]** Helper function to pull out average expression and
fraction of cells for genes of interest.

## Usage

``` r
rs_compute_cluster_expr_stats(f_path_gene, gene_indices, clusters)
```

## Arguments

- f_path_gene:

  Path to the `counts_genes.bin` file.

- gene_indices:

  Integer vector. 0-indexed(!) positions of the genes to include.

- clusters:

  List. A list that contains within the cell indices of the clusters of
  interest (0-indexed).

## Value

A list with two matrices

- mean - The average expression. Shape genes x clusters.

- frac - The fraction of cells expressing. Shape genes x clusters.
