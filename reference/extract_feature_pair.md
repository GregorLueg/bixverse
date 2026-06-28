# Extract a pair of features for scatter / hex plots

Extracts two features into a wide data.table with `feature_1` and
`feature_2` value columns, ready for a scatter or hex plot. Each feature
may carry a `_rna` or `_adt` suffix to choose its modality independently
(e.g. `"ENSG00000167286_rna"` against `"CD3_adt"`); features without a
suffix fall back to `modality`. For `SingleCells` / `MetaCells` only RNA
exists, so an `_adt` feature there errors via
[`extract_gene_expression()`](https://gregorlueg.github.io/bixverse/reference/extract_gene_expression.md).

## Usage

``` r
extract_feature_pair(
  object,
  feature_1,
  feature_2,
  obs_cols = NULL,
  scale = FALSE,
  clip = NULL,
  modality = c("rna", "adt")
)
```

## Arguments

- object:

  A single cell class.

- feature_1:

  String. First feature, optionally `_rna` / `_adt` suffixed.

- feature_2:

  String. Second feature, optionally `_rna` / `_adt` suffixed.

- obs_cols:

  Optional character vector. Obs columns to attach (e.g. to colour the
  scatter).

- scale:

  Boolean. Whether to z-score the expression values per feature.

- clip:

  Optional numeric. Clip z-scores if `scale = TRUE`.

- modality:

  String. Fallback modality for unsuffixed features. One of
  `c("rna", "adt")`.

## Value

A data.table with `cell_id`, `feature_1`, `feature_2` and any requested
obs columns. The original feature labels are stored in a `features`
attribute as `c(feature_1, feature_2)`.
