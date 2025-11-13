# Calculate the percentage of gene sets in the cells

This function allows to calculate for example the proportion of
mitochondrial genes, or ribosomal genes in the cells for QC purposes.

## Usage

``` r
rs_sc_get_gene_set_perc(
  f_path_cell,
  gene_set_idx,
  cell_indices,
  streaming,
  verbose
)
```

## Arguments

- f_path_cell:

  String. Path to the `counts_cells.bin` file.

- gene_set_idx:

  Named list with integer(!) positions (0-indexed!) as elements of the
  genes of interest.

- cell_indices:

  Integer. The indices of the cells for which to calculate the
  proportions. (0-indexed!)

- streaming:

  Boolean. Shall the data be worked on in chunks.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with the percentages of counts per gene set group detected in the
cells.
