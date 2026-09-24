# Calculates the gene statistics for a set of cell groups and genes

**\[experimental\]** Helper function to extract data for dot plots
and/or heatmaps.

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

  Integer vector. Gene index positions to return (0-indexed!).

- group_ids:

  Integer vector. Group of each cell in `cell_indices`, as an index into
  `group_levels` (0-indexed!). Same length as `cell_indices`.

- group_levels:

  Character vector. The group labels.

## Value

A list with the following elements:

- grp_label - The group labels, i.e. `group_levels`.

- mean_exp - Mean normalised expression per gene and group over all
  cells of the group (zeros included), row-major (genes x groups).

- perc_exp - Fraction (`[0, 1]`) of cells in the group with a non-zero
  count, row-major (genes x groups).
