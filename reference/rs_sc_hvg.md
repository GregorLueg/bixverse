# Calculate the percentage of gene sets in the cells

**\[experimental\]** This function identifies highly variable genes with
the three methods known in Seurat.

## Usage

``` r
rs_sc_hvg(
  f_path_gene,
  hvg_method,
  cell_indices,
  loess_span,
  binning,
  n_bins,
  clip_max,
  streaming,
  verbose
)
```

## Arguments

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- hvg_method:

  String. Which HVG detection method to use. One of
  `c("vst", "meanvarbin", "dispersion")`.

- cell_indices:

  Integer positions (0-indexed!) that defines the cells to keep.

- loess_span:

  Numeric. The span parameter for the loess function.

- binning:

  String. The binning strategy for the `meanvarbin` method. One of
  `c("equal_width", "equal_frequency")`.

- n_bins:

  Integer. Number of bins for the `meanvarbin` method.

- clip_max:

  Optional clipping number. Defaults to `sqrt(no_cells)` if not
  provided.

- streaming:

  Boolean. Shall the genes be streamed in to reduce memory pressure.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the highly variable genes. If `hvg_method == "vst"`, the
following elements can be found:

- mean - The average expression of the gene.

- var - The variance of the gene.

- var_exp - The expected variance of the gene.

- var_std - The standardised variance of the gene.

For the other two methods, these elements can be found:

- mean - The average expression of the gene.

- dispersion - The dispersion of the gene

- dispersion_scaled - The scaled dispersion per bin per gene

- bin - The bin of the gene
