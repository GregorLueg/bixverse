# Helper to extract single cell counts for several genes

Helper to extract single cell counts for several genes

## Usage

``` r
rs_extract_several_genes_plots(f_path, cell_indices, gene_indices, scale, clip)
```

## Arguments

- f_path:

  String. Path to the `counts_genes.bin` file.

- cell_indices:

  Integer positions (0-indexed!) that defines the cells to keep.

- gene_indices:

  Integer. Gene index position to return (0-indexed!).

- scale:

  Boolean. Shall the normalised counts be scaled.

- clip:

  Optional float. Clipping for the Z-scores if scale is set to `TRUE`

## Value

A list of dense vectors of the normalised counts.
