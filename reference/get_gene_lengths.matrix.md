# Get gene set lengths for a matrix

Get gene set lengths for a matrix

## Usage

``` r
# S3 method for class 'matrix'
get_gene_lengths(x, species = c("human", "mouse", "rat"), ...)
```

## Arguments

- x:

  Numerical matrix. Assumes genes x samples as format and Ensembl
  identifiers as gene ids.

- species:

  String. One of `c("human", "mouse", "rat")`.

- ...:

  Additional parameters. Not in use atm.

## Value

Named numeric representing the gene lengths.

## Examples

``` r
if (FALSE) { # \dontrun{
# the matrix method, dispatched on Ensembl rownames
counts <- matrix(
  1:4,
  nrow = 2,
  dimnames = list(c("ENSG00000141510", "ENSG00000012048"), c("s1", "s2"))
)
get_gene_lengths.matrix(counts, species = "human")
} # }
```
