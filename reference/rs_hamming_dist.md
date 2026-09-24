# Calculates the Hamming distance between categorical columns

**\[experimental\]** Calculates the pairwise Hamming distance between
the columns, i.e. the fraction of rows in which two columns differ.

## Usage

``` r
rs_hamming_dist(x)
```

## Arguments

- x:

  Integer matrix. The integers represent the factor data.

## Value

The Hamming distance matrix between the columns, values in `[0, 1]`.
