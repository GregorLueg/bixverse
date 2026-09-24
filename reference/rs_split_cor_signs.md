# Helper function to split correlation matrices by sign

**\[experimental\]** All `1` if every off-diagonal correlation is
non-negative, all `-1` if every one is non-positive, otherwise a
graph-based split into two groups.

## Usage

``` r
rs_split_cor_signs(data)
```

## Arguments

- data:

  The correlation matrix to split by sign.

## Value

An integer vector of 1 and -1, one per column of `data`.
