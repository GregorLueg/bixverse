# Calculate the effect size

Calculate the effect size

## Usage

``` r
hedges_g_dge(
  meta_data,
  main_contrast,
  normalised_counts,
  contrast_list = NULL,
  .verbose = TRUE
)
```

## Arguments

- meta_data:

  data.table. The meta information about the experiment in which the
  contrast info can be found.

- main_contrast:

  String. Which column contains the main groups you want to calculate
  the Hedge's G effect for. Every permutation of the groups will be
  tested if `contrast_list` is `NULL`.

- normalised_counts:

  Numeric Matrix. The normalised count matrix.

- contrast_list:

  String vector or NULL. Optional string vector of contrast formatted as
  `"contrast1-contrast2"`. Default NULL will create all contrasts
  automatically.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

A data.table with the effect sizes and standard errors based on the
Hedge's G effect size for the groups.

## Examples

``` r
# Hedge's G on log CPM counts for the case vs control contrast
syn <- synthetic_bulk_cor_matrix()
meta <- data.table::data.table(
  sample_id = colnames(syn$counts),
  case_control = rep(c("case", "control"), each = 50)
)
norm_counts <- edgeR::cpm(syn$counts, log = TRUE)
res <- hedges_g_dge(
  meta_data = meta,
  main_contrast = "case_control",
  normalised_counts = norm_counts,
  .verbose = FALSE
)
head(res)
#>    effect_sizes standard_errors gene_id     combination
#>           <num>           <num>  <char>          <char>
#> 1:   -0.5265904       0.2034367  gene_1 case_vs_control
#> 2:   -0.2279491       0.2006485  gene_2 case_vs_control
#> 3:   -0.4553987       0.2025758  gene_3 case_vs_control
#> 4:   -0.4393655       0.2023986  gene_4 case_vs_control
#> 5:   -0.3488589       0.2015155  gene_5 case_vs_control
#> 6:   -0.4848591       0.2029173  gene_6 case_vs_control
```
