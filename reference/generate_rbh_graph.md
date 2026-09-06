# Generate an RBH graph.

This function will generate an RBH graph based on set similarity between
gene modules. You have the option to use an overlap coefficient instead
of Jaccard similarity and to specify a minimum similarity.

## Usage

``` r
generate_rbh_graph(
  object,
  minimum_similarity,
  k_best = 1L,
  overlap_coefficient = FALSE,
  spearman = FALSE
)
```

## Arguments

- object:

  The underlying class, see
  [`RbhGraph()`](https://gregorlueg.github.io/bixverse/reference/RbhGraph.md).

- minimum_similarity:

  The minimum similarity to create an edge.

- k_best:

  Integer. Number of best neighbours to consider. If set to `1L`, this
  behaves as the traditional reciprocal best hit. If you set this to
  `3L` you consider edges if the modules is in the top 3 best modules by
  similarity for each other.

- overlap_coefficient:

  Boolean. Shall the overlap coefficient be used instead of Jaccard
  similarity. Only relevant if the underlying class is set to set
  similarity.

- spearman:

  Boolean. Shall Spearman correlation be used. Only relevant if the
  underlying class is set to correlation-based similarity.

## Value

The class with added properties.

## Examples

``` r
# Jaccard-based reciprocal best hits between two module sets
set.seed(123)
modules <- data.table::data.table(
  origin = rep(c("set_a", "set_b"), each = 20),
  module = rep(c("m1", "m2", "m3", "m4"), each = 10),
  gene = unlist(replicate(4, sample(letters, 10), simplify = FALSE))
)
object <- RbhGraph(
  modules,
  rbh_type = "set",
  dataset_col = "origin",
  module_col = "module",
  value_col = "gene"
)
object <- generate_rbh_graph(object, minimum_similarity = 0)
head(get_rbh_res(object))
#>    origin target origin_modules target_modules similiarity combined_origin
#>    <char> <char>         <char>         <char>       <num>          <char>
#> 1:  set_a  set_b             m1             m3   0.2500000        set_a_m1
#> 2:  set_a  set_b             m2             m4   0.3333333        set_a_m2
#>    combined_target
#>             <char>
#> 1:        set_b_m3
#> 2:        set_b_m4
```
