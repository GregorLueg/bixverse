# Bayesian PCA on a matrix with missing values

**\[experimental\]** Port of `bpca()` from pcaMethods. Variational Bayes
with an ARD prior per component, so superfluous components shrink
towards zero. Deterministic, the start comes from an SVD.

## Usage

``` r
rs_bpca(x, bpca_params, verbose)
```

## Arguments

- x:

  Numeric matrix. Rows = samples, columns = features. `NA` marks a
  missing value.

- bpca_params:

  List. The BPCA parameters, see
  [`params_bpca()`](https://gregorlueg.github.io/bixverse/reference/params_bpca.md).

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with:

- scores - Samples x `n_pcs`.

- loadings - Features x `n_pcs`, not orthonormal.

- r2_cum - Cumulative R^2 per component on the observed entries.

- centre - Column centres that were subtracted.

- scale - Column scales that were divided out.

- completed - `x` with the missing entries imputed.

- noise_var - Residual variance, `1 / tau`.

- n_iter - Variational steps run.

- converged - Whether `tol` was reached before `max_iter`.

## References

Oba, et al., Bioinformatics, 2003; Stacklies, et al., Bioinformatics,
2007
