# Pseudo-bulk a set of cells (sparse)

This function will return a sparse matrix of
`length(cell_indices_ls) x number of genes` (in list form in CSR). The
function has the option to return the sum of the sum of the raw counts
or the average of the normalised counts.

## Usage

``` r
rs_pseudobulk_cells_sparse(f_path, cell_indices_ls, assay, verbose)
```

## Arguments

- f_path:

  String. Path to the `counts_cells.bin` file.

- cell_indices_ls:

  List. Must contains 0-indexed positions of the cells to aggregate per
  element.

- assay:

  String. One of `c("raw", "norm")`. Which counts to normalise.

- verbose:

  Controls verbosity of the function

## Value

A list with the following elements (easy to convert into CSR in R)

- indptr - The index pointers (representing cells)

- indices - The indices (representing genes)

- data - The pseudo-bulked data

- nrow - Number of rows (i.e., pseudo-bulked samples)

- ncol - Number of columns
