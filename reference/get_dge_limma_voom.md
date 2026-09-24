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
object <- calculate_dge_limma(
  object,
  contrast_column = "case_control",
  .verbose = FALSE
)
head(get_dge_limma_voom(object))
#>     gene_id      logFC        CI.L        CI.R   AveExpr         t     P.Value
#>      <char>      <num>       <num>       <num>     <num>     <num>       <num>
#> 1: gene_589  0.2246372  0.08408076  0.36519360 10.469973  3.169891 0.002011843
#> 2:  gene_25 -0.6156898 -1.03578573 -0.19559385  9.230661 -2.906878 0.004475873
#> 3: gene_894  0.2106366  0.05376385  0.36750934 10.165859  2.663175 0.008992482
#> 4: gene_691  0.1905472  0.04817933  0.33291502 10.089603  2.654630 0.009208069
#> 5: gene_649 -0.2152072 -0.37767245 -0.05274193 10.033452 -2.627300 0.009929367
#> 6:  gene_24 -0.2383683 -0.42102215 -0.05571450 10.247703 -2.588410 0.011043868
#>    adj.P.Val         B        contrast subgroup
#>        <num>     <num>          <char>   <lgcl>
#> 1: 0.7339682 -2.136281 case_vs_control       NA
#> 2: 0.7339682 -3.252944 case_vs_control       NA
#> 3: 0.7339682 -3.072213 case_vs_control       NA
#> 4: 0.7339682 -3.113507 case_vs_control       NA
#> 5: 0.7339682 -3.175555 case_vs_control       NA
#> 6: 0.7339682 -3.150449 case_vs_control       NA
```
