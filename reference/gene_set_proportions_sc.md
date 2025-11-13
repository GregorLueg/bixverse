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
  streaming = FALSE,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- gene_set_list:

  A named list with each element containing the gene identifiers of that
  set. These should be the same as `get_gene_names(object)`!

- streaming:

  Boolean. Shall the cells be streamed in. Useful for larger data sets
  where you wish to avoid loading in the whole data. Default to `FALSE`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

It will add the columns based on the names in the `gene_set_list` to the
obs table.
