# Generate pseudo-bulked matrices

This is a wrapper function to generate pseudo-bulked counts from a list
of cell identifiers. You have the option to return sparse or dense
matrices and if you wish to sum raw counts or calculate the mean of the
normalised counts.

## Usage

``` r
get_pseudobulked_sc(
  object,
  cell_list,
  return_format = c("dense", "sparse"),
  assay = c("raw", "norm"),
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- cell_list:

  Named list.

- return_format:

  String. One of `c("dense", "sparse")`.

- assay:

  String. One of `c("raw", "norm")`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

Pending on your setting in return_format a dense matrix or sparse CSR
matrix with aggregated cells x genes.
