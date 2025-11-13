# Get gene set lengths for a bulk_coexp class

Get gene lengths for a bulk_coexp object by extracting the count matrix
and delegating to the matrix method.

## Arguments

- x:

  An object of class `bulk_coexp`.

- species:

  String. One of `c("human", "mouse", "rat")`.

- ...:

  Additional parameters. Not in use atm.

## Value

Named numeric representing the gene lengths.
