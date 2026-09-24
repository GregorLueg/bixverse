# Generates synthetic single cell counts with a planted ambient profile

**\[experimental\]** Builds the fixture CellSweep is tested against.
Every real barcode is a two-component multinomial: a planted fraction
`alpha` of its library comes from the soup, the rest from its own cell
type profile. Empty droplets are pure soup at a much smaller library
size, which is what the ambient profile is estimated off. Real barcodes
come first in the matrix, empty droplets after, and cell types are
assigned round-robin over the real barcodes.

The soup is cell type one plus flat background rather than a mixture of
every profile. A soup sitting in the span of the cell type profiles
makes the contamination fraction unidentifiable.

## Usage

``` r
rs_synthetic_sc_cellsweep_data(
  n_real,
  n_empty,
  n_genes,
  n_celltypes,
  n_markers,
  marker_weight,
  ambient_dominance,
  alpha_mean,
  alpha_sd,
  real_lib_size,
  empty_lib_size,
  seed
)
```

## Arguments

- n_real:

  Integer. Number of real barcodes.

- n_empty:

  Integer. Number of empty droplets. At least 30, and the ambient
  profile gets noisy well above that.

- n_genes:

  Integer. Number of genes.

- n_celltypes:

  Integer. Number of cell types.

- n_markers:

  Integer. Width of each cell type's marker block. Blocks are disjoint
  and laid out from the first gene.

- marker_weight:

  Float. Enrichment of a marker gene over background in its own cell
  type's profile. Must exceed 1.

- ambient_dominance:

  Float. Fraction of the soup coming from the first cell type. The
  remainder is flat background.

- alpha_mean:

  Float. Mean planted ambient fraction across real barcodes.

- alpha_sd:

  Float. Spread of the planted ambient fraction.

- real_lib_size:

  Integer. Expected library size of a real barcode.

- empty_lib_size:

  Integer. Expected library size of an empty droplet.

- seed:

  Integer. For reproducibility.

## Value

A list with the following items.

- data - Integer vector. Non-zero counts of the CSR matrix.

- indptr - Integer vector. Row pointers of the CSR matrix.

- indices - Integer vector. 0-indexed(!) gene positions.

- nrow - Integer. Number of barcodes, real plus empty.

- ncol - Integer. Number of genes.

- cell_type_indices - Integer vector. 0-indexed(!) cell type per real
  barcode. Empty droplets have none.

- is_empty - Logical vector over all barcodes.

- alpha_true - Numeric vector. Planted ambient fraction per real
  barcode.

- ambient_true - Numeric vector. The soup, summing to one.

- celltype_profiles_true - Numeric vector. Cell type profiles, row-major
  `n_celltypes x n_genes`, each row summing to one.
