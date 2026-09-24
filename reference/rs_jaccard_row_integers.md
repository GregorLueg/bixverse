# Calculate rapidly Jaccard similarities between rows

**\[experimental\]** Helper function to quickly calculate the Jaccard
similarity between matching rows of the two matrices. Each row is
treated as a set of integers (duplicates removed).

## Usage

``` r
rs_jaccard_row_integers(data_1, data_2)
```

## Arguments

- data_1:

  Integer matrix. The first matrix to compare.

- data_2:

  Integer matrix. The second matrix to compare. Needs the same number of
  rows as `data_1`.

## Value

The Jaccard similarity averaged over the rows.
