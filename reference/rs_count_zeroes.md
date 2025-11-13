# Helper to get zero stats from a given matrix

Calculates in a single matrix pass the total number of zeroes, the row
zeroes and column zeroes.

## Usage

``` r
rs_count_zeroes(x)
```

## Arguments

- x:

  Numeric matrix. The matrix for which to calculate the total zeroes,
  column and row zeroes.

## Value

A list with:

- total_zeroes - Total zeroes in the matrix.

- row_zeroes - Vector with number of zeroes per row.

- col_zeroes - Vector with number of zeroes per column.
