# Wrapper function for the branch cell selection parameters

Parameters controlling which cells
[`run_gene_trends_sc()`](https://gregorlueg.github.io/bixverse/reference/run_gene_trends_sc.md)
assigns to each branch. The threshold on a fate's probability is an
expanding quantile over the pseudotime-sorted cells, made monotone with
a cumulative maximum, so a fate's bar can only rise as differentiation
proceeds. The defaults are the reference ones.

## Usage

``` r
params_sc_branch_selection(q = 0.01, eps = 0.01, resolution = 500L)
```

## Arguments

- q:

  Numeric. Upper-tail quantile of the fate probability used as the
  threshold. Defaults to `0.01`.

- eps:

  Numeric. Slack subtracted from the threshold before the comparison.
  Defaults to `0.01`.

- resolution:

  Integer. Number of pseudotime buckets, capped at the cell count.
  Defaults to `500L`.

## Value

A named flat list with all branch selection parameters.

## References

Setty, et al., Nat. Biotechnol., 2019.
