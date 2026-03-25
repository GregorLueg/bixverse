# Calculates the gene statistics for a set of cell groups and genes

Helper function to extract data for dot plots and/or heatmaps.

## Usage

``` r
rs_extract_grouped_gene_stats(
  f_path,
  cell_indices,
  gene_indices,
  group_ids,
  group_levels
)
```

## Arguments

- f_path:

  String. Path to the `counts_genes.bin` file.

- cell_indices:

  Integer positions (0-indexed!) that defines the cells to keep.

- gene_indices:

  Integer. Gene index position to return (0-indexed!).

- group_ids:

  Integer. The levels of the data. (0-indexed!)

- group_levels:

  String. Name of the factors.

## Value

A list with the following elements:

- grp_label - The label of that group

- mean_exp - Vector of mean expression values in row major (genes x
  n_levels)

- perc_exp - Vector of proportions of cells with expression in row major
  (genes x n_levels)
