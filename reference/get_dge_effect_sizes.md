# Return the effect size results

Getter function to extract the Effect size results from the
[`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md)
class.

## Usage

``` r
get_dge_effect_sizes(object)
```

## Arguments

- object:

  `BulkDge` class.

## Value

Returns the effect size results. (If found.)

## Examples

``` r
# Hedge's G effect sizes per contrast
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
object <- calculate_dge_hedges(
  object,
  contrast_column = "case_control",
  .verbose = FALSE
)
head(get_dge_effect_sizes(object))
#>    effect_sizes standard_errors gene_id     combination subgroup
#>           <num>           <num>  <char>          <char>   <lgcl>
#> 1:   -0.4541874       0.2046602  gene_1 case_vs_control       NA
#> 2:   -0.2674739       0.2029738  gene_2 case_vs_control       NA
#> 3:   -0.4256828       0.2043474  gene_3 case_vs_control       NA
#> 4:   -0.3793411       0.2038811  gene_4 case_vs_control       NA
#> 5:   -0.2601789       0.2029254  gene_5 case_vs_control       NA
#> 6:   -0.4551123       0.2046707  gene_6 case_vs_control       NA
```
