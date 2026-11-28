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
the reference, or name one group to only get that arm. Either way,
everything is calculated in a single pass over the genes.

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
  method = NULL,
  alternative = c("greater", "less", "twosided"),
  min_prop = 0.05,
  modality = c("rna", "atac"),
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

  Optional string. One of `c("wilcox", "lr")`: the Wilcoxon rank sum
  test, or Seurat's logistic regression test (`LRDETest()`) with the
  library size as the latent variable. `NULL` picks `"lr"` on ATAC data
  and `"wilcox"` otherwise.

- alternative:

  String. Test alternative. One of `c("twosided", "greater", "less")`.
  This function will default to `"greater"`, i.e., genes upregulated in
  the reference group.

- min_prop:

  Numeric. The minimum proportion of cells that need to express the gene
  in at least one of the groups. Applied once, globally, so every
  comparison's FDR is calculated over the same gene set.

- modality:

  String. One of `c("rna", "atac")`. `"atac"` tests the features of the
  ATAC store of a `SingleCellsMultiModal`; on a `SingleCellsAtac` both
  name its own store. Results are then keyed by `feature_id`.

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
#> 1: cell_type_1 gene_01 0.9880239    0.9531267 0.9520958  0.9531267 0.9541575
#> 2: cell_type_1 gene_02 0.9880239    0.9434446 0.9398312  0.9434446 0.9470580
#> 3: cell_type_1 gene_03 1.0000000    0.9699717 0.9675348  0.9699717 0.9724085
#> 4: cell_type_1 gene_04 0.9580838    0.9179307 0.9086916  0.9179307 0.9271697
#> 5: cell_type_1 gene_05 0.8383234    0.8323234 0.8302915  0.8323234 0.8343554
#> 6: cell_type_1 gene_06 0.9760479    0.9415449 0.9400656  0.9415449 0.9430242
#>    worst_rival min_rank      simes_p    simes_fdr        max_p    max_p_fdr
#>         <char>    <int>        <num>        <num>        <num>        <num>
#> 1: cell_type_3        4 6.581099e-47 8.226374e-46 1.432128e-46 1.790160e-45
#> 2: cell_type_3        5 1.615884e-45 1.615884e-44 3.255928e-44 2.713274e-43
#> 3: cell_type_3        2 1.378589e-50 3.446472e-49 1.178490e-49 2.946226e-48
#> 4: cell_type_2        7 1.081619e-41 7.725847e-41 1.125127e-38 8.036625e-38
#> 5: cell_type_2        9 2.159089e-27 1.199494e-26 5.643505e-27 3.135280e-26
#> 6: cell_type_3        5 9.874756e-45 8.228963e-44 2.109971e-44 2.109971e-43

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
