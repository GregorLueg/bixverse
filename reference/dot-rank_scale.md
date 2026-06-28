# Rank-based min-max scaling

Returns `rank(x) / max(rank(x))` with average ties and NAs sent to the
bottom. Used to map differential-expression statistics into `[0, 1]`
before aggregation.

## Usage

``` r
.rank_scale(x)
```

## Arguments

- x:

  Numeric vector.

## Value

Numeric vector in `[0, 1]`.
