# Generates synthetic data with a planted multicellular programme

**\[experimental\]** Builds the fixture DIALOGUE is tested against.
Every cell type gets its own noise and its own sample-level nuisance
factors; only the first feature column and the planted genes carry the
shared per-sample latent, so anything found beyond that is spurious.
Cells are laid out contiguously by cell type and, within a cell type, by
sample.

The counts are a scaled copy of the normalised layer rather than a draw
from a count model. That keeps the planted signal clean, which is the
point of a fixture.

## Usage

``` r
rs_synthetic_sc_dialogue_data(
  n_samples,
  cells_per_sample,
  n_cell_types,
  n_features,
  n_sample_features,
  n_genes,
  n_planted,
  seed
)
```

## Arguments

- n_samples:

  Integer. Samples the experiment spans. DIALOGUE needs at least 5.

- cells_per_sample:

  Integer. Cells per sample per cell type.

- n_cell_types:

  Integer. Number of cell types. Must be at least 2.

- n_features:

  Integer. Feature columns per cell type. Must be at least 2.

- n_sample_features:

  Integer. Feature columns carrying a per-sample component. The first of
  those is the shared programme, the rest are cell-type-specific
  nuisance; anything past this count is pure noise.

- n_genes:

  Integer. Number of genes.

- n_planted:

  Integer. Planted genes per cell type. Cell type `t` owns genes
  `t * n_planted` to `(t + 1) * n_planted - 1` (0-indexed), so the
  blocks have to fit into `n_genes`.

- seed:

  Integer. Random seed for reproducibility.

## Value

A list with the following items.

- data - The synthetic raw counts, CSR over cells.

- indptr - The index pointers of the cells.

- indices - The gene indices for the given cells.

- nrow - Number of cells.

- ncol - Number of genes.

- cell_type_indices - List of integer vectors. 0-indexed(!) global cell
  positions per cell type.

- features - List of numeric matrices, one per cell type, rows aligned
  to `cell_type_indices`.

- sample_ids - Integer vector. 0-indexed(!) sample code per cell.

- quality - Numeric vector. Quality covariate per cell. Pure noise.

- latent - Numeric vector. The per-sample latent the planted programme
  follows.

- planted - List of integer vectors. 0-indexed(!) planted gene positions
  per cell type.

## References

Jerby-Arnon & Regev, Nature Biotechnology, 2022
