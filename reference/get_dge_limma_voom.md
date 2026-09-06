# Return the Limma Voom results

Getter function to extract the Limma Voom results from the
[`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md)
class.

## Usage

``` r
get_dge_limma_voom(object)
```

## Arguments

- object:

  `BulkDge` class.

## Value

Returns the Limma Voom results. (If found.)

## Examples

``` r
# topTable results for every contrast that was fitted
syn <- synthetic_bulk_cor_matrix()
meta <- data.table::data.table(
  sample_id = colnames(syn$counts),
  case_control = rep(c("case", "control"), each = 50)
)
object <- BulkDge(raw_counts = syn$counts, meta_data = meta)
object <- qc_bulk_dge(object, group_col = "case_control", .verbose = FALSE)
object <- normalise_bulk_dge(
  object,
  group_col = "case_control",
  .verbose = FALSE
)
#> calcNormFactors has been renamed to normLibSizes
object <- calculate_dge_limma(
  object,
  contrast_column = "case_control",
  .verbose = FALSE
)
head(get_dge_limma_voom(object))
#>     gene_id      logFC       CI.L        CI.R   AveExpr         t     P.Value
#>      <char>      <num>      <num>       <num>     <num>     <num>       <num>
#> 1:  gene_25 -0.6397658 -1.0498106 -0.22972104  9.230661 -3.094264 0.002538899
#> 2:  gene_24 -0.2766326 -0.4614834 -0.09178174 10.247703 -2.967901 0.003726954
#> 3: gene_649 -0.2592920 -0.4354871 -0.08309685 10.033452 -2.918521 0.004317106
#> 4: gene_754 -0.2463135 -0.4141150 -0.07851200 10.039387 -2.911119 0.004412630
#> 5: gene_461 -0.2253523 -0.3791305 -0.07157408 10.468820 -2.906262 0.004476363
#> 6:  gene_30 -0.3027633 -0.5170655 -0.08846100  9.576105 -2.801844 0.006068301
#>    adj.P.Val         B        contrast subgroup
#>        <num>     <num>          <char>   <lgcl>
#> 1: 0.5821864 -3.118852 case_vs_control       NA
#> 2: 0.5821864 -2.709253 case_vs_control       NA
#> 3: 0.5821864 -2.876284 case_vs_control       NA
#> 4: 0.5821864 -2.883107 case_vs_control       NA
#> 5: 0.5821864 -2.702562 case_vs_control       NA
#> 6: 0.5821864 -3.231806 case_vs_control       NA
```
