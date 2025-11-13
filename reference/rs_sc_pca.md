# Calculates PCA for single cell

Helper function that will calculate the PCA for the specified highly
variable genes. Has the option to use randomised SVD for faster solving
of the PCA.

## Usage

``` r
rs_sc_pca(
  f_path_gene,
  no_pcs,
  random_svd,
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

- return_scaled:

  Boolean. Shall the scaled data be returned.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with with the following items

- scores - The samples projected on the PCA space.

- loadings - The loadings of the features for the PCA.

- scaled - The scaled matrix if you set return_scaled to `TRUE`.
