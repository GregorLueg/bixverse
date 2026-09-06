# Wrapper for a Limma Voom analysis

Wrapper function to run Limma Voom workflows.

## Usage

``` r
run_limma_voom(
  meta_data,
  main_contrast,
  dge_list,
  contrast_list = NULL,
  co_variates = NULL,
  quantile_norm = FALSE,
  ...,
  .verbose = TRUE
)
```

## Arguments

- meta_data:

  data.table. The meta information about the experiment in which the
  contrast info (and potential co-variates) can be found.

- main_contrast:

  String. Which column contains the main groups you want to test
  differential gene expression with the Limma-Voom workflow for.

- dge_list:

  DGEList, see
  [`edgeR::DGEList()`](https://rdrr.io/pkg/edgeR/man/DGEList.html).

- contrast_list:

  String vector or NULL. Optional string vector of contrast formatted as
  `"contrast1-contrast2"`. Default NULL will create all contrasts
  automatically.

- co_variates:

  String or NULL. Optional co-variates you wish to consider during model
  fitting.

- quantile_norm:

  Boolean. Shall the counts be also quantile-normalised. Defaults to
  `FALSE`.

- ...:

  Additional parameters to forward to
  [`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html) or
  [`limma::voom()`](https://rdrr.io/pkg/limma/man/voom.html).

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A data.table with all the DGE results from
[`limma::topTable()`](https://rdrr.io/pkg/limma/man/toptable.html) for
the identified contrast pairs.

## Examples

``` r
# voom fit and topTable results for the single case vs control contrast
syn <- synthetic_bulk_cor_matrix()
meta <- data.table::data.table(
  sample_id = colnames(syn$counts),
  case_control = rep(c("case", "control"), each = 50)
)
dge_list <- edgeR::normLibSizes(edgeR::DGEList(counts = syn$counts))
res <- run_limma_voom(
  meta_data = meta,
  main_contrast = "case_control",
  dge_list = dge_list,
  .verbose = FALSE
)
head(res)
#>     gene_id      logFC        CI.L        CI.R   AveExpr         t     P.Value
#>      <char>      <num>       <num>       <num>     <num>     <num>       <num>
#> 1: gene_589  0.2299554  0.09099159  0.36891915 10.476457  3.281322 0.001405009
#> 2:  gene_25 -0.6558752 -1.07387792 -0.23787240  9.202591 -3.111350 0.002401557
#> 3: gene_691  0.2078402  0.06467766  0.35100274 10.103653  2.878770 0.004842709
#> 4: gene_894  0.2222151  0.06532984  0.37910034 10.175698  2.808654 0.005938488
#> 5: gene_612  0.1948064  0.05683337  0.33277942 10.436486  2.799728 0.006093182
#> 6:  gene_69 -0.4358784 -0.77662943 -0.09512729  8.571075 -2.536500 0.012675264
#>    adj.P.Val         B        contrast
#>        <num>     <num>          <char>
#> 1:   0.74139 -1.920476 case_vs_control
#> 2:   0.74139 -3.035099 case_vs_control
#> 3:   0.74139 -2.781126 case_vs_control
#> 4:   0.74139 -2.850416 case_vs_control
#> 5:   0.74139 -2.753920 case_vs_control
#> 6:   0.74139 -3.841576 case_vs_control
```
