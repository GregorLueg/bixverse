# Probabilistic PCA on a matrix with missing values

**\[experimental\]** Port of `ppca()` from pcaMethods. Fits the
principal subspace by EM on the observed entries only, then imputes the
missing entries from it. Use it to get a PCA or a completed matrix out
of data with sporadic `NA`s, for example proteomics or metabolomics
intensities.

## Usage

``` r
run_ppca(x, ppca_params = params_ppca(), .verbose = TRUE)
```

## Arguments

- x:

  Numeric matrix. Rows = samples, columns = features. `NA` (or `NaN`)
  marks a missing value. No row or column may be entirely missing and
  observed values must be finite.

- ppca_params:

  List, see
  [`params_ppca()`](https://gregorlueg.github.io/bixverse/reference/params_ppca.md).
  `n_pcs` must not exceed `min(dim(x))`.

- .verbose:

  Boolean or integer. Verbosity.

## Value

A list with:

- scores - Samples x `n_pcs`.

- loadings - Features x `n_pcs`, orthonormal.

- r2_cum - Cumulative R^2 per component on the completed matrix.

- centre - Column centres that were subtracted.

- scale - Column scales that were divided out.

- completed - `x` with the missing entries imputed. Observed entries are
  untouched.

- noise_var - Residual variance outside the subspace.

- n_iter - EM iterations run.

- converged - Whether `tol` was reached before `max_iter`.

## Details

The initial loadings are drawn from `ppca_params$seed` with the Rust
RNG, so results agree with pcaMethods at convergence, not iterate by
iterate. PPCA loadings are orthonormal. Signs of the components are
arbitrary.

## References

Roweis, NIPS, 1998; Tipping and Bishop, J R Stat Soc B, 1999; Stacklies,
et al., Bioinformatics, 2007

## Examples

``` r
set.seed(42L)
x <- matrix(rnorm(50L * 10L), nrow = 50L, ncol = 10L)
x[sample(length(x), 50L)] <- NA
ppca_res <- run_ppca(x, .verbose = FALSE)
head(ppca_res$scores)
#>              PC1        PC2
#> [1,] -2.23444601 -2.0828126
#> [2,]  0.10485233  0.7263406
#> [3,]  0.04822818 -0.6237992
#> [4,] -0.13365587  1.3769826
#> [5,] -0.20892364 -0.4253578
#> [6,] -0.01437090 -0.4680165
```
