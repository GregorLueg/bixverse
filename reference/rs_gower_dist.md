# Calculates the Gower distance for a given matrix

**\[experimental\]** Calculates the pairwise Gower distance between the
rows. Continuous columns contribute the range-normalised absolute
difference, categorical columns a simple mismatch.

## Usage

``` r
rs_gower_dist(x, is_cat)
```

## Arguments

- x:

  Numerical matrix. Converted matrix of continuous and categorical
  variables as numerical values. Rows = samples, columns = features.

- is_cat:

  Logical vector of length `ncol(x)`. Which of the columns represent
  categorical values.

## Value

The Gower distance matrix between the rows, values in `[0, 1]`.
