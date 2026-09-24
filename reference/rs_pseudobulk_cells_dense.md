# Pseudo-bulk a set of cells (dense)

**\[experimental\]** This function will return a dense matrix of
`length(cell_indices_ls) x number of genes`. The function has the option
to return the sum of the raw counts or the average of the normalised
counts.

## Usage

``` r
rs_pseudobulk_cells_dense(f_path, cell_indices_ls, assay, verbose)
```

## Arguments

- f_path:

  String. Path to the `counts_cells.bin` file.

- cell_indices_ls:

  List. Each element contains the 0-indexed positions of the cells to
  aggregate.

- assay:

  String. One of `c("raw", "norm")`. `"raw"` sums the raw counts,
  `"norm"` averages the normalised counts. Unrecognised values fall back
  to `"raw"`.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A dense numerical matrix of pseudo-bulked samples x genes.
