# Wrapper function for gene trend parameters

Parameters controlling the landmark Gaussian process that
[`run_gene_trends_sc()`](https://gregorlueg.github.io/bixverse/reference/run_gene_trends_sc.md)
fits per branch. The kernel is a Matern-5/2 one and the prediction grid
doubles as the landmark set.

The defaults come from the reference and are prior-dominated. Palantir's
pseudotime is min-max scaled to `[0, 1]`, so a `length_scale` of `1.0`
spans the entire domain and a `sigma` of `1.0` sits at roughly the
signal scale of log-normalised expression. The posterior will flatten
genuine transient structure and resolve almost any gene into a smooth
monotone or single-peaked curve. That is a presentation choice, not
inference. Shorten `length_scale` before believing a bump.

## Usage

``` r
params_sc_gene_trends(
  resolution = 500L,
  weighting = c("hard_mask", "fate_probability"),
  length_scale = 1,
  sigma = 1,
  jitter = 1e-06,
  max_jitter_retries = 3L,
  chunk_size = 2048L
)
```

## Arguments

- resolution:

  Integer. Grid points per branch. Kept at the default even when a
  branch holds fewer cells, as the reference does. Defaults to `500L`.

- weighting:

  String. One of `c("hard_mask", "fate_probability")`. With
  `"hard_mask"` every selected cell enters its branch's fit with equal
  weight, which is what the reference does. With `"fate_probability"`
  every cell enters every fit weighted by its fate probability, which is
  more defensible: a cell at 0.6 is not a member. Defaults to
  `"hard_mask"`.

- length_scale:

  Numeric. Matern-5/2 length scale. Defaults to `1.0`.

- sigma:

  Numeric. Noise standard deviation. Defaults to `1.0`.

- jitter:

  Numeric. Added to the landmark covariance diagonal before the
  Cholesky. Defaults to `1e-6`.

- max_jitter_retries:

  Integer. Times the jitter is raised and the Cholesky retried before
  giving up. Defaults to `3L`.

- chunk_size:

  Integer. Training points held at once when accumulating the
  cross-covariance. Defaults to `2048L`.

## Value

A named flat list with all gene trend parameters. The Gaussian process
hyperparameters sit at the same level as the trend ones, not in a nested
block.

## References

Setty, et al., Nat. Biotechnol., 2019.
