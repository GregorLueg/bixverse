# Pseudo-bulk a set of cells (dense)

This function will return a dense matrix of
`length(cell_indices_ls) x number of genes`. The function has the option
to return the sum of the sum of the raw counts or the average of the
normalised counts.

## Usage

``` r
rs_pseudobulk_cells_dense(f_path, cell_indices_ls, assay, verbose)
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

A dense matrix with the pseudo-bulked data.
