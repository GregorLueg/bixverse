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

## Examples

``` r
# markers that hold up against every rival cell type
sc <- demo_single_cells()
res <- find_specific_markers_sc(
  sc,
  column_of_interest = "cell_grp",
  .verbose = FALSE
)
head(res$summary)
#>        ref_grp gene_id  prop_ref median_auroc min_auroc mean_auroc max_auroc
#>         <char>  <char>     <num>        <num>     <num>      <num>     <num>
#> 1: cell_type_1 gene_01 0.9880239    0.9531087 0.9520777  0.9531087 0.9541396
#> 2: cell_type_1 gene_02 0.9880239    0.9434354 0.9397771  0.9434354 0.9470938
#> 3: cell_type_1 gene_03 1.0000000    0.9699715 0.9674987  0.9699715 0.9724444
#> 4: cell_type_1 gene_04 0.9580838    0.9179306 0.9087095  0.9179306 0.9271517
#> 5: cell_type_1 gene_05 0.8383234    0.8323145 0.8302556  0.8323145 0.8343734
#> 6: cell_type_1 gene_06 0.9760479    0.9415538 0.9400476  0.9415538 0.9430600
#>    worst_rival min_rank      simes_p    simes_fdr        max_p    max_p_fdr
#>         <char>    <int>        <num>        <num>        <num>        <num>
#> 1: cell_type_3        4 6.632426e-47 8.290532e-46 1.443307e-46 1.804134e-45
#> 2: cell_type_3        5 1.589201e-45 1.589201e-44 3.333359e-44 2.777799e-43
#> 3: cell_type_3        2 1.354582e-50 3.386456e-49 1.197853e-49 2.994632e-48
#> 4: cell_type_2        7 1.089790e-41 7.784218e-41 1.116568e-38 7.975489e-38
#> 5: cell_type_2        9 2.145260e-27 1.191811e-26 5.714230e-27 3.174572e-26
#> 6: cell_type_3        5 9.715047e-45 8.095872e-44 2.126360e-44 2.126360e-43

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
