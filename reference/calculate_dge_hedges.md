# Calculates the Hedge's G effect size

This function will calculate the Hedge's G effect size on the normalised
counts. Should batch-corrected counts be found, these will be used. At a
minimum you will need to provide `contrast_column` that can be found in
the meta-data. If you do not provide a vector of contrasts that you wish
to test for, every permutation of groups represented in that column will
be tested against each other.

## Usage

``` r
calculate_dge_hedges(
  object,
  contrast_column,
  contrast_list = NULL,
  filter_column = NULL,
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

- .verbose:

  Controls verbosity of the function.

## Value

Returns the class with additional data added to the outputs.

## Examples

``` r
# effect sizes over every contrast in the case_control column
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
