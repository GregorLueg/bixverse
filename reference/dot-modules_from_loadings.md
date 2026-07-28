# Derive sparse module membership from a loading matrix

Turns a `gene x k` loading matrix from a matrix factorisation (ICA, NMF,
DGRDL) into a membership table by keeping the tails of each component's
loading distribution.

The important property is that membership is **not** exclusive: a gene
that loads strongly on three components appears three times. That is the
whole point of a factorisation, and it is why an argmax assignment is
the wrong tool. Genes that fail the threshold on every component appear
not at all, giving a real background category rather than forcing every
feature into some module.

Two thresholding rules, selected via `membership_params`:

- `"zscore"` - robust standardisation per component, centred on the
  median and scaled by the MAD, keeping `abs(z) > cutoff`. No
  distributional assumption beyond rough symmetry.

- `"fdr"` - two-sided p-values against a Normal null fitted per
  component by median and MAD, Benjamini-Hochberg adjusted, keeping
  `padj < fdr`.

## Usage

``` r
.modules_from_loadings(
  loadings,
  membership_params = params_module_membership()
)
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
