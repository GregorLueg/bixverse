# Resolve cell and gene selection for NMF on SingleCells

Resolve cell and gene selection for NMF on SingleCells

## Usage

``` r
.resolve_sc_nmf_selection(object, cell_ids, gene_ids)
```

## Arguments

- object:

  `SingleCells` class.

- cell_ids:

  Optional string vector. The cells to include.

- gene_ids:

  Optional string vector. The genes to include.

## Value

A list with the following items:

- cell_indices - The Rust-based cell indices.

- gene_indices - The Rust-based gene indices.

- cell_ids - The cell identifiers.

- gene_ids - The gene identifiers.
