# Calculate HVG per batch

Batch-aware highly variable gene detection. Calculates HVG statistics
separately for each batch, allowing for downstream selection strategies
such as union of top genes per batch.

## Usage

``` r
rs_sc_hvg_batch_aware(
  f_path_gene,
  hvg_method,
  cell_indices,
  batch_labels,
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

- batch_labels:

  Integer vector (0-indexed!) defining batch membership for each cell.
  Must be same length as `cell_indices`.

- loess_span:

  Numeric. The span parameter for the loess function.

- binning:

  String. The binning strategy for the `meanvarbin` method. One of
  `c("equal_width", "equal_frequency")`.

- n_bins:

  Integer. Number of bins for the `meanvarbin` method.

- clip_max:

  Optional clipping number. Defaults to `sqrt(no_cells)` per batch if
  not provided.

- streaming:

  Boolean. Shall the genes be streamed in to reduce memory pressure.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with HVG statistics concatenated across all batches. For
`hvg_method == 'vst'`, the following elements can be found:

- mean - The average expression of each gene in each batch.

- var - The variance of each gene in each batch.

- var_exp - The expected variance of each gene in each batch.

- var_std - The standardised variance of each gene in each batch.

- batch - Batch index for each gene (length = n_genes \* n_batches).

- gene_idx - Gene index for each entry (0-indexed, length = n_genes \*
  n_batches).

For the other methods

- mean - The average expression of each gene in each batch.

- dispersion - The dispersion of the gene in each batch.

- dispersion_scaled - The scaled dispersion per bin per gene in each
  batch.

- bin - The bin of the gene in each batch.

- batch - Batch index for each gene (length = n_genes \* n_batches).

- gene_idx - Gene index for each entry (0-indexed, length = n_genes \*
  n_batches).
