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
  streaming = FALSE,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- top_n_vals:

  Integer. The Top N thresholds to test.

- streaming:

  Boolean. Shall the cells be streamed in. Useful for larger data sets
  where you wish to avoid loading in the whole data. Default to `FALSE`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

It will add the columns based on the names in the `gene_set_list` to the
obs table.
