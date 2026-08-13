# Find markers that are specific to a cell group

This function scores one group of cells against every other group of a
column separately, instead of against all of them pooled together. That
difference matters for markers. A pooled test is dominated by whichever
rival contributes the most cells, so a gene that is high in the
reference and just as high in one small rival still comes out looking
like a clean marker. Here the gene has to hold up against every rival,
and the per-gene summaries (`min_auroc`, `median_auroc`, `min_rank`)
tell you whether it does.

Leave `reference_group` as `NULL` to run every group of the column as
the reference in turn, or name one group to only get that arm.

The summaries rank on AUROC rather than the p-value on purpose. Group
sizes vary a lot in practice and p-values scale with the group sizes, so
a large rival would otherwise crowd out a small one regardless of effect
size.

## Usage

``` r
find_specific_markers_sc(
  object,
  column_of_interest,
  reference_group = NULL,
  method = "wilcox",
  alternative = c("greater", "less", "twosided"),
  min_prop = 0.05,
  downsampling = TRUE,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- column_of_interest:

  String. The column in the obs table holding the groups, e.g. a cell
  type annotation or a clustering.

- reference_group:

  Optional string. The group to use as the reference. If `NULL`, every
  group of `column_of_interest` is used as the reference in turn.

- method:

  String. Which method to use for the calculations of the DGE. At the
  moment the only option is `"wilcox"`, but the parameter is reserved
  for future features.

- alternative:

  String. Test alternative. One of `c("twosided", "greater", "less")`.
  This function will default to `"greater"`, i.e., genes upregulated in
  the reference group.

- min_prop:

  Numeric. The minimum proportion of cells that need to express the gene
  in at least one of the groups. Applied once, globally, so every
  comparison's FDR is calculated over the same gene set.

- downsampling:

  Boolean. If any group exceeds 100,000 cells, a random subsample of
  100,000 cells is used for it. The subsample is drawn once per group,
  so a group is represented by the same cells in every arm.

- seed:

  Integer. Seed that is used for the downsampling.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

A `ScSpecificMarkers` class with the following elements

- summary - data.table. Per gene and reference group, the summaries
  across all rivals: `prop_ref`, `median_auroc`, `min_auroc`,
  `mean_auroc`, `max_auroc`, `worst_rival` (the rival achieving
  `min_auroc`), `min_rank` (best AUROC rank the gene reaches against any
  single rival), the Simes-combined p-value with its FDR, and the
  maximum p-value with its FDR.

- per_comparison - data.table. Per gene, reference group and rival, the
  underlying `auroc`, `lfc`, `prop_ref`, `prop_rival`, `z_scores`,
  `p_values` and `fdr`.

- params - List. The parameters the run used.

## References

Soneson and Robinson, Nat Methods, 2018; Lun, et al., F1000Research,
2016
