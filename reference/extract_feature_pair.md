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
  modality = c("rna", "adt"),
  layer = c("norm", "magic")
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

- layer:

  String. One of `c("norm", "magic")`, forwarded to
  [`extract_gene_expression()`](https://gregorlueg.github.io/bixverse/reference/extract_gene_expression.md).
  Applies to both features, and an `_adt` suffixed one will error under
  `"magic"`.

## Value

A data.table with `cell_id`, `feature_1`, `feature_2` and any requested
obs columns. The original feature labels are stored in a `features`
attribute as `c(feature_1, feature_2)`.

## Examples

``` r
# two genes side by side, ready for a scatter
sc <- demo_single_cells(prepped = FALSE)
genes <- get_gene_names(sc)[1:2]
head(extract_feature_pair(sc, genes[1], genes[2], obs_cols = "cell_grp"))
#> Key: <cell_id>
#>     cell_id feature_1 feature_2    cell_grp
#>      <char>     <num>     <num>      <char>
#> 1: cell_001  6.222656  6.761719 cell_type_1
#> 2: cell_002  4.113281  4.511719 cell_type_2
#> 3: cell_003  4.300781  4.300781 cell_type_3
#> 4: cell_004  5.757812  5.757812 cell_type_1
#> 5: cell_005  0.000000  0.000000 cell_type_2
#> 6: cell_006  0.000000  4.234375 cell_type_3

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
