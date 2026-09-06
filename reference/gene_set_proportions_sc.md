# Calculate the proportions of reads for specific gene sets

This is a helper function that calculates proportions of reads belonging
to given gene sets. This can be used for example for the calculation of
percentage mitochondrial reads per cell. These will be automatically
added to the obs table

## Usage

``` r
gene_set_proportions_sc(
  object,
  gene_set_list,
  streaming = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- gene_set_list:

  A named list with each element containing the gene identifiers of that
  set. These should be the same as `get_gene_names(object)`!

- streaming:

  Optional Boolean. Shall the data be streamed in. Useful for larger
  data sets where you wish to avoid loading in the whole data. If
  `NULL`, will automatically detect.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

It will add the columns based on the names in the `gene_set_list` to the
obs table.

## Examples

``` r
# read proportion of a gene set, the mitochondrial percentage pattern
sc <- demo_single_cells(prepped = FALSE)
sc <- gene_set_proportions_sc(
  sc,
  gene_set_list = list(set_a = c("gene_01", "gene_02", "gene_03")),
  .verbose = FALSE
)
head(unlist(sc[["set_a"]]))
#>      set_a1      set_a2      set_a3      set_a4      set_a5      set_a6 
#> 0.161870509 0.021021022 0.016949153 0.283667624 0.000000000 0.009070295 

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
