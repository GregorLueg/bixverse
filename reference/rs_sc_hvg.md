# Calculate the percentage of gene sets in the cells

This function allows to calculate for example the proportion of
mitochondrial genes, or ribosomal genes in the cells for QC purposes.

## Usage

``` r
rs_sc_hvg(
  f_path_gene,
  hvg_method,
  cell_indices,
  loess_span,
  clip_max,
  streaming,
  verbose
)
```

## Arguments

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- hvg_method:

  String. Which HVG detection method to use. Options are
  `c("vst", "meanvarbin", "dispersion")`. So far, only the first is
  implemented.

- cell_indices:

  Integer positions (0-indexed!) that defines the cells to keep.

- loess_span:

  Numeric. The span parameter for the loess function.

- clip_max:

  Optional clipping number. Defaults to `sqrt(no_cells)` if not
  provided.

- streaming:

  Boolean. Shall the genes be streamed in to reduce memory pressure.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with the percentages of counts per gene set group detected in the
cells:

- mean - The average expression of the gene.

- var - The variance of the gene.

- var_exp - The expected variance of the gene.

- var_std - The standardised variance of the gene.
