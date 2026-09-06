# Helper function to transform the Rust-exported sparse matrices into R ones

Helper function to transform the Rust-exported sparse matrices into R
ones

## Usage

``` r
sparse_list_to_mat(ls)
```

## Arguments

- ls:

  List. Needs to represent the (column) sparse data.

## Value

The sparseMatrix from the data.

## Examples

``` r
# 3 x 3 identity in the Rust-side CSC representation
csc <- list(
  data = c(1, 1, 1),
  indices = c(0L, 1L, 2L),
  indptr = c(0L, 1L, 2L, 3L),
  nrow = 3L,
  ncol = 3L,
  cs_type = "csc"
)
sparse_list_to_mat(csc)
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>           
#> [1,] 1 . .
#> [2,] . 1 .
#> [3,] . . 1
```
