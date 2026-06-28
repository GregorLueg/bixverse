# Calculates PCA for single cell

**\[experimental\]** Helper function that will calculate the PCA for the
specified highly variable genes. You have the option to do mean
centering, variance normalisation and/or apply the new proposed
transformation `PFlogPF` from Booeshaghi, et al.

## Usage

``` r
rs_sc_pca(
  f_path_gene,
  f_path_cell,
  no_pcs,
  pca_params,
  cell_indices,
  gene_indices,
  seed,
  return_scaled,
  verbose
)
```

## Arguments

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- f_path_cell:

  String. Path to the `counts_cells.bin` file. Used if you wish to use
  the PFlogPF transformation.

- no_pcs:

  Integer. Number of PCs to calculate.

- pca_params:

  Named list. Contains the parameters to use for this PCA run.

- cell_indices:

  Integer. The cell indices to use. (0-indexed!)

- gene_indices:

  Integer. The gene indices to use. (0-indexed!)

- seed:

  Integer. Random seed for the randomised SVD.

- return_scaled:

  Boolean. Shall the scaled data be returned.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with with the following items

- scores - The samples projected on the PCA space.

- loadings - The loadings of the features for the PCA.

- singular_values - The singular values for the PCA.

- scaled - The scaled matrix if you set return_scaled to `TRUE`.

## References

Booeshaghi, et al., bioRxive, 2026.
