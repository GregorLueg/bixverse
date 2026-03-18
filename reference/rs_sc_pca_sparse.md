# Calculates sparse PCA for single cell

Helper function that will calculate sparse PCA without scaling the data.
This has the advantage that you avoid creating a large dense matrix due
to scaling; however, it has the disadvantage that the first PC will be
heavily influenced by average expression. If random_svd is set to
`FALSE`, Lanczos iterations will be used to solve the SVD; if random_svd
is set to `TRUE`, the randomised version will be used with
multiplication of the initial sparse matrix with a much smaller random
dense matrix, avoiding holding a large dense matrix in memory.

## Usage

``` r
rs_sc_pca_sparse(
  f_path_gene,
  no_pcs,
  random_svd,
  cell_indices,
  gene_indices,
  seed,
  verbose
)
```

## Arguments

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- no_pcs:

  Integer. Number of PCs to calculate.

- random_svd:

  Boolean. Shall randomised SVD be used.

- cell_indices:

  Integer. The cell indices to use. (0-indexed!)

- gene_indices:

  Integer. The gene indices to use. (0-indexed!)

- seed:

  Integer. Random seed for the randomised SVD.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with with the following items

- scores - The samples projected on the PCA space (solved via sparse
  SVD).

- loadings - The loadings of the features for the PCA (solved via sparse
  SVD).

- singular_values - The singular values for the PCA (solved via sparse
  SVD).
