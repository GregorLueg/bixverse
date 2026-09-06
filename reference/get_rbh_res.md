# Get the RBH results

Pulls out the RBH results if you ran
[`generate_rbh_graph()`](https://gregorlueg.github.io/bixverse/reference/generate_rbh_graph.md)

## Usage

``` r
get_rbh_res(object)
```

## Arguments

- object:

  The underlying class
  [`RbhGraph()`](https://gregorlueg.github.io/bixverse/reference/RbhGraph.md).

## Value

The data.table with the RBH result if found, otherwise NULL.

## Examples

``` r
# reciprocal best hits between modules of two data sets
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
