# Derive sparse module membership from a loading matrix

Turns a `gene x k` loading matrix from a matrix factorisation (ICA, NMF,
DGRDL) into a membership table by keeping the tails of each component's
loading distribution.

Two thresholding rules, selected via `membership_params`:

- `"zscore"` - standardisation per component, keeping `abs(z) > cutoff`.
  No distributional assumption beyond rough symmetry.

- `"fdr"` - two-sided p-values against a Normal null fitted per
  component, Benjamini-Hochberg adjusted, keeping `padj < fdr`.

The standardisation itself is controlled by `membership_params$scaling`:
`"robust"` centres and scales by the median and MAD, `"standard"` by the
mean and standard deviation. The latter is stricter and keeps fewer
genes on skewed loadings, which is common for NMF.

## Usage

``` r
modules_from_loadings(loadings, membership_params = params_module_membership())
```

## Arguments

- loadings:

  Numeric matrix. `gene x k`, with row and column names.

- membership_params:

  List. See
  [`params_module_membership()`](https://gregorlueg.github.io/bixverse/reference/params_module_membership.md).

## Value

A data.table with columns `gene`, `module_id`, `loading`, `sign` and the
per-component score (`z` or `padj` depending on the method). One row per
surviving (gene, component) pair, ordered by component then by
descending absolute loading.
