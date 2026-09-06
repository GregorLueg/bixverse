# Get the gene lengths

Get the gene lengths

## Usage

``` r
get_gene_lengths(x, species = c("human", "mouse", "rat"), ...)

## S7 method for class <bixverse::BulkCoExp>
get_gene_lengths(x, species = c("human", "mouse", "rat"), ...)
```

## Arguments

- x:

  Object to extract gene lengths from. Can be a matrix or BulkCoExp
  object.

- species:

  String. One of `c("human", "mouse", "rat")`.

- ...:

  Additional parameters passed to methods.

## Examples

``` r
if (FALSE) { # \dontrun{
# median transcript length per Ensembl gene, queried from Ensembl
counts <- matrix(
  1:4,
  nrow = 2,
  dimnames = list(c("ENSG00000141510", "ENSG00000012048"), c("s1", "s2"))
)
get_gene_lengths(counts, species = "human")
} # }
```
