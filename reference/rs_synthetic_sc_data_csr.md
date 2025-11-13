# Generate synthetic single cell data (h5ad type)

This function generates pseudo data to test single cell functions in
form of the h5ad version, with cells = rows and genes = columns. The
data is encoded in CSR.

## Usage

``` r
rs_synthetic_sc_data_csr(n_genes, n_cells, min_genes, max_genes, max_exp, seed)
```

## Arguments

- n_genes:

  Integer. Number of genes you wish to have in the synthetic data.

- n_cells:

  Integer. Number of cells you wish to have in the synthetic data.

- min_genes:

  Integer. Minimum number of genes expressed per cell.

- max_genes:

  Integer. Maximum number of genes expressed per cell.

- max_exp:

  Upper bound in terms of expression. Expression values will be sampled
  from `1:max_exp`.

- seed:

  Integer. Seed for reproducibility purposes.

## Value

The list with the synthetic data with the following items:

- data - The synthetic counts

- row_ptrs - The row pointers

- col_indices - The column indices

- nrow - Number of rows (cells)

- ncol - Number of cols (genes)
