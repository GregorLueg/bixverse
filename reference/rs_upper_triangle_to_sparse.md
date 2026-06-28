# Generate sparse data from an upper triangle

**\[experimental\]** This function takes the values from an upper
triangle matrix the shift and the nrows/ncols and returns a list.

## Usage

``` r
rs_upper_triangle_to_sparse(value, shift, n, cs_type)
```

## Arguments

- value:

  Numeric vector. The upper triangle values.

- shift:

  Boolean. Was the matrix shifted up (false = diagonal included; true
  diagonal not incldued).

- n:

  Integer. The number of columns/rows in the symmetric matrix.

- cs_type:

  String. One of `c("csr", "csc")`. Which type of list to return.

## Value

A list containing:

- data - A vector of lists with the elements. (Related to the way Robj
  are stored in Rust.)

- row_indices - A vector of integers with the row indices.

- col_ptr - A vector of integers with the column pointers.
