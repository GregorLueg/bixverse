# Calculate the critical value

**\[experimental\]** Same as
[`rs_critval()`](https://gregorlueg.github.io/bixverse/reference/rs_critval.md),
but the values are taken from the upper triangle (diagonal excluded) of
a symmetric matrix.

## Usage

``` r
rs_critval_mat(mat, iters, alpha, seed)
```

## Arguments

- mat:

  Numeric matrix. The symmetric matrix with all of the values.

- iters:

  Integer. Size of the bootstrap sample.

- alpha:

  Float. The alpha value. For example, 0.001 would return the value
  exceeded by roughly 0.1 percent of the bootstrap sample.

- seed:

  Integer. For reproducibility purposes

## Value

The critical value for the given parameters.
