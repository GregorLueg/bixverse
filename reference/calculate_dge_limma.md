# Calculates the Limma Voom DGE

This function will apply the Limma Voom DGE workflow. At a minimum you
will need to provide `contrast_column` that can be found in the
meta-data. If you do not provide a vector of contrasts that you wish to
test for, every permutation of groups represented in that column will be
tested against each other.

## Usage

``` r
calculate_dge_limma(
  object,
  contrast_column,
  contrast_list = NULL,
  filter_column = NULL,
  co_variates = NULL,
  limma_params = params_limma_voom(),
  .verbose = TRUE
)
```

## Arguments

- object:

  The underlying class, see
  [`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md).

- contrast_column:

  String. The contrast column in which the groupings are stored. Needs
  to be found in the meta_data within the properties.

- contrast_list:

  Optional string vector. A vectors that contains the contrast formatted
  as `"contrast1-contrast2"`. Default `NULL` will create all possible
  contrast automatically.

- filter_column:

  Optional String. If there is a column you wish to use as sub
  groupings, this can be provided here. An example could be different
  sampled tissues and you wish to run the DGE analyses within each
  tissue separately in the data.

- co_variates:

  Optional string vector. Any co-variates you wish to consider during
  the Limma Voom modelling.

- limma_params:

  List. The limma parameters, see
  [`params_limma_voom()`](https://gregorlueg.github.io/bixverse/reference/params_limma_voom.md).

- .verbose:

  Controls verbosity of the function.

## Value

Returns the class with additional data added to the outputs.

## Examples

``` r
# limma voom over every contrast in the case_control column
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
