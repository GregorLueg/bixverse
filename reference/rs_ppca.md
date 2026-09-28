# Probabilistic PCA on a matrix with missing values

**\[experimental\]** Port of `ppca()` from pcaMethods. Fits the
principal subspace by EM on the observed entries only and imputes the
missing ones from it. The start is drawn from `seed` with the Rust RNG,
so results match pcaMethods at convergence, not iterate by iterate.

## Usage

``` r
rs_ppca(x, ppca_params, verbose)
```

## Arguments

- x:

  Numeric matrix. Rows = samples, columns = features. `NA` marks a
  missing value.

- ppca_params:

  List. The PPCA parameters, see
  [`params_ppca()`](https://gregorlueg.github.io/bixverse/reference/params_ppca.md).

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with:

- scores - Samples x `n_pcs`.

- loadings - Features x `n_pcs`, orthonormal.

- r2_cum - Cumulative R^2 per component on the completed matrix.

- centre - Column centres that were subtracted.

- scale - Column scales that were divided out.

- completed - `x` with the missing entries imputed.

- noise_var - Residual variance outside the subspace.

- n_iter - EM iterations run.

- converged - Whether `tol` was reached before `max_iter`.

## References

Roweis, NIPS, 1998; Tipping and Bishop, J R Stat Soc B, 1999; Stacklies,
et al., Bioinformatics, 2007
