# Calculate the proportions of reads for the Top N genes

This is a helper function that calculates proportions of reads to the
Top N genes by expression in a given cell. High values here can indicate
low complexity, quality cells. The values will be automatically added to
the obs table.

## Usage

``` r
top_genes_perc_sc(
  object,
  top_n_vals = c(25L, 50L, 100L),
  streaming = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- top_n_vals:

  Integer. The Top N thresholds to test.

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
# share of a cell's reads taken by its top 5 and top 10 genes
sc <- demo_single_cells(prepped = FALSE)
sc <- top_genes_perc_sc(sc, top_n_vals = c(5L, 10L), .verbose = FALSE)
head(unlist(sc[["top_5_genes_percentage"]]))
#> top_5_genes_percentage1 top_5_genes_percentage2 top_5_genes_percentage3 
#>               0.5863310               0.6096096               0.6924940 
#> top_5_genes_percentage4 top_5_genes_percentage5 top_5_genes_percentage6 
#>               0.7106017               0.6787879               0.7664399 

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
