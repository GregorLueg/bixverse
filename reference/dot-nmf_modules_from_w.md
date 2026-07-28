# Derive gene-to-module assignments from an NMF loadings matrix

Each gene is assigned to the component whose absolute loading is
highest. Thin wrapper over
[`.modules_from_loadings()`](https://gregorlueg.github.io/bixverse/reference/dot-modules_from_loadings.md).
NMF loadings are non-negative, so the `"auto"` tail setting resolves to
an upper-tail test.

## Usage

``` r
.nmf_modules_from_w(w, membership_params = params_module_membership())
```

## Arguments

- w:

  Numeric matrix. Gene loadings (features x k).

- membership_params:

  List. See
  [`params_module_membership()`](https://gregorlueg.github.io/bixverse/reference/params_module_membership.md).

## Value

A data.table with one row per surviving (gene, component) pair.
