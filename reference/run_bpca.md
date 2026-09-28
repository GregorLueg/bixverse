# Bayesian PCA on a matrix with missing values

**\[experimental\]** Port of `bpca()` from pcaMethods. Variational Bayes
with an ARD prior per component, so components that the data do not
support shrink towards zero rather than fitting noise. Missing entries
are imputed from the fit.

## Usage

``` r
run_bpca(x, bpca_params = params_bpca(), .verbose = TRUE)
```

## Arguments

- x:

  Numeric matrix. Rows = samples, columns = features. `NA` (or `NaN`)
  marks a missing value. No row or column may be entirely missing and
  observed values must be finite.

- bpca_params:

  List, see
  [`params_bpca()`](https://gregorlueg.github.io/bixverse/reference/params_bpca.md).
  `n_pcs` must not exceed `min(dim(x))`.

- .verbose:

  Boolean or integer. Verbosity.

## Value

A list with:

- scores - Samples x `n_pcs`.

- loadings - Features x `n_pcs`, not orthonormal.

- r2_cum - Cumulative R^2 per component on the observed entries.

- centre - Column centres that were subtracted.

- scale - Column scales that were divided out.

- completed - `x` with the missing entries imputed. Observed entries are
  untouched.

- noise_var - Residual variance, `1 / tau`.

- n_iter - Variational steps run.

- converged - Whether `tol` was reached before `max_iter`.

## Details

Deterministic, the start comes from an SVD. BPCA loadings are not
orthonormal, and `r2_cum` is computed on the observed entries only, as
in pcaMethods. Signs of the components are arbitrary.

## References

Oba, et al., Bioinformatics, 2003; Stacklies, et al., Bioinformatics,
2007

## Examples

``` r
set.seed(42L)
x <- matrix(rnorm(50L * 10L), nrow = 50L, ncol = 10L)
x[sample(length(x), 50L)] <- NA
bpca_res <- run_bpca(x, .verbose = FALSE)
head(bpca_res$scores)
#>                PC1           PC2
#> [1,]  1.564589e-13  2.147490e-14
#> [2,] -2.172977e-14 -1.312883e-14
#> [3,]  1.171312e-14  1.760569e-14
#> [4,] -3.228109e-14 -2.846022e-14
#> [5,]  1.998706e-14  1.640825e-15
#> [6,]  1.146295e-14  1.106884e-14
```
