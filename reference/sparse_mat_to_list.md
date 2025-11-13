# Helper function to transform R sparse matrices to list

Helper function to transform R sparse matrices to list

## Usage

``` r
sparse_mat_to_list(sparse_mat)
```

## Arguments

- sparse_mat:

  Sparse matrix. The matrix to transform into a list.

## Value

A list with the following elements

- indptr - Index pointers of the sparse data.

- indices - Indices of the data.

- data - The underlying data.

- format - String that defines if the data is CSR or CSC.

- nrow - The number of rows.

- ncol - The number of columns.
