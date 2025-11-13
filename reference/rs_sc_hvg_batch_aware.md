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
  clip_max,
  streaming,
  verbose
)
```

## Arguments

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- hvg_method:

  String. Which HVG detection method to use. Currently only `"vst"` is
  implemented for batch-aware mode.

- cell_indices:

  Integer positions (0-indexed!) that defines the cells to keep.

- batch_labels:

  Integer vector (0-indexed!) defining batch membership for each cell.
  Must be same length as `cell_indices`.

- loess_span:

  Numeric. The span parameter for the loess function.

- clip_max:

  Optional clipping number. Defaults to `sqrt(no_cells)` per batch if
  not provided.

- streaming:

  Boolean. Shall the genes be streamed in to reduce memory pressure.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with HVG statistics concatenated across all batches:

- mean - The average expression of each gene in each batch.

- var - The variance of each gene in each batch.

- var_exp - The expected variance of each gene in each batch.

- var_std - The standardised variance of each gene in each batch.

- batch - Batch index for each gene (length = n_genes \* n_batches).

- gene_idx - Gene index for each entry (0-indexed, length = n_genes \*
  n_batches).
