# Generate sparse data from an upper triangle

**\[experimental\]** This function takes the values from an upper
triangle matrix, the shift and the nrows/ncols and returns the full
symmetric matrix as a compressed sparse list.

## Usage

``` r
rs_upper_triangle_to_sparse(value, shift, n, cs_type)
```

## Arguments

- value:

  Numeric vector. The upper triangle values.

- shift:

  Boolean. Was the matrix shifted up (`FALSE` = diagonal included;
  `TRUE` = diagonal not included).

- n:

  Integer. The number of columns/rows in the symmetric matrix.

- cs_type:

  String. One of `c("csr", "csc")`. Which type of list to return. Other
  values raise an error.

## Value

A list containing:

- data - Numeric vector with the non-zero values.

- indptr - Integer vector with the index pointers.

- indices - Integer vector with the 0-based indices.

- nrow - Number of rows.

- ncol - Number of columns.

- cs_type - `"csr"` or `"csc"`.
