# Resolve what a NEBULA Wald test is testing

The single cell counterpart of
[`.resolve_tested()`](https://gregorlueg.github.io/bixverse/reference/dot-resolve_tested.md).
NEBULA tests one coefficient or one contrast, not a set of them, so this
returns a scalar `coef` or a single weight vector rather than edgeR's
`n_contrasts` form.

## Usage

``` r
.resolve_tested_sc(design, coef = NULL, contrast = NULL)
```

## Arguments

- design:

  Numeric matrix. The design matrix of cells x coefficients.

- coef:

  Optional integer or character. The coefficient to report, as a 1-based
  column position or a column name.

- contrast:

  Optional numeric vector. One weight per column of `design`.

## Value

A named list holding either `coef` (0-based) or `contrast`.
