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
