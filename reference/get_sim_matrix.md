# Get the similarity matrix

Get the similarity matrix

## Usage

``` r
get_sim_matrix(object, as_data_table = FALSE, .verbose = TRUE)
```

## Arguments

- object:

  `OntologySim class`. See
  [`OntologySim()`](https://gregorlueg.github.io/bixverse/reference/OntologySim.md).

- as_data_table:

  Boolean. Shall the data be returned as a long data.table.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

Returns the semantic similarity data.table from the class

## Examples

``` r
# the Wang similarity matrix back out of the class
onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f"),
  type = c("part_of", "part_of", "part_of", "is_a", "is_a")
)
onto_obj <- OntologySim(onto, .verbose = FALSE)
onto_obj <- calculate_wang_sim_onto(
  onto_obj,
  weights = c(part_of = 0.8, is_a = 0.6),
  .verbose = FALSE
)
round(get_sim_matrix(onto_obj, .verbose = FALSE), 3)
#>       f     d     b     a     c     e
#> f 1.000 0.470 0.625 0.400 0.796 0.428
#> d 0.470 1.000 0.764 0.477 0.590 0.558
#> b 0.625 0.764 1.000 0.643 0.764 0.742
#> a 0.400 0.477 0.643 1.000 0.477 0.481
#> c 0.796 0.590 0.764 0.477 1.000 0.558
#> e 0.428 0.558 0.742 0.481 0.558 1.000
```
