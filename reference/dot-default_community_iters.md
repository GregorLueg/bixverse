# Default number of community detection iterations

`5L` for Leiden and `2L` for Louvain, which gave the best modularity
across synthetic graphs.

## Usage

``` r
.default_community_iters(algorithm)
```

## Arguments

- algorithm:

  String. One of `c("leiden", "louvain")`.

## Value

Integer.
