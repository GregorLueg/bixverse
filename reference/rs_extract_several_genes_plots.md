# Helper to extract single cell counts for several genes

**\[experimental\]** Extract the normalised single cell counts of
several genes at once.

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

  Integer vector. Gene index positions to return (0-indexed!).

- scale:

  Boolean. Shall the normalised counts be z-scored per gene across the
  selected cells.

- clip:

  Optional float. Clips the Z-scores to `[-clip, clip]`. Only used if
  `scale = TRUE`.

## Value

A list of numerical vectors, one per gene in `gene_indices`, each with
one normalised value per cell in `cell_indices`.
