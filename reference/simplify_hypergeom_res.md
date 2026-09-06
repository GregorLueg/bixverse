# Simplify gene set results via ontologies

This function provides an interface to simplify overenrichment results
based on ontological information of the pathway origin (typical use case
is to simplify gene ontology results). To do so, the function will
calculate the Wang similarity and keep within a set of highly similar
terms the one with the lowest fdr. Should there be terms with the same
fdr, the function will keep the most specific term within the ontology.

## Usage

``` r
simplify_hypergeom_res(res, parent_child_dt, weights, min_sim = 0.7)
```

## Arguments

- res:

  data.table. The enrichment results. Needs to have the columns
  `c("gene_set_name", "fdr")`.

- parent_child_dt:

  data.table. The data.table with column parent and child. You also need
  to have a type column for the Wang similarity to provide the weights
  for the relationships.

- weights:

  Named numeric. The relationship of type to weight for this specific
  edge. For example `c("part_of" = 0.8, "is_a" = 0.6)`.

- min_sim:

  Float between 0 and 1. The minimum similarity that the terms need to
  have.

## Value

data.table with enrichment results.

## Examples

``` r
# collapse redundant terms of a toy ontology by Wang similarity
onto <- data.table::data.table(
  parent = c("a", "b", "b", "c"),
  child = c("b", "c", "d", "e"),
  type = c("is_a", "is_a", "part_of", "is_a")
)
res <- data.table::data.table(
  gene_set_name = c("b", "c", "d", "e"),
  fdr = c(0.001, 0.001, 0.4, 0.02)
)
simplify_hypergeom_res(
  res,
  onto,
  weights = c("is_a" = 0.8, "part_of" = 0.6),
  min_sim = 0.6
)
#>    gene_set_name   fdr
#>           <char> <num>
#> 1:             e  0.02
```
