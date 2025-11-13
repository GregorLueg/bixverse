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
