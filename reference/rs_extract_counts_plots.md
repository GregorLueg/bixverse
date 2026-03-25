# Helper to extract single cell counts as a dense vector for plotting

Helper to extract single cell counts as a dense vector for plotting

## Usage

``` r
rs_extract_counts_plots(f_path, cell_indices, gene_index, norm, scale, clip)
```

## Arguments

- f_path:

  String. Path to the `counts_genes.bin` file.

- cell_indices:

  Integer positions (0-indexed!) that defines the cells to keep.

- gene_index:

  Integer. Gene index position to return (0-indexed!).

- norm:

  Boolean. Shall normalised counts be returned.

- scale:

  Boolean. Shall the normalised counts be scaled.

- clip:

  Optional float. Clipping for the Z-scores if scale is set to `TRUE`

## Value

The dense vector of expression values for this gene.
