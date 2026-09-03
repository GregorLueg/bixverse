# Resolve what an edgeR-style test is testing

Turns the user-facing `coef` / `contrast` pair into the elements the
Rust side reads off the parameter list. `coef` accepts 1-based column
positions or column names of the design and is shifted to the 0-based
indices Rust wants. `contrast` is flattened column-major with its column
count alongside.

Supplying neither tests the last column of the design, which is what
edgeR and limma default to.

## Usage

``` r
.resolve_tested(design, coef = NULL, contrast = NULL)
```

## Arguments

- design:

  Numeric matrix. The design matrix of samples x coefficients.

- coef:

  Optional integer or character. The coefficient(s) to drop from the
  null model.

- contrast:

  Optional numeric vector or matrix. Weights over the coefficients, one
  entry (or row) per column of `design`.

## Value

A named list holding either `coef` or `contrast` plus `n_contrasts`.
