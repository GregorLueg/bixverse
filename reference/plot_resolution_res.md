# Plot the resolution results.

Plots the resolution results (if they can be found in the class). The
x-axis reflects the different resolutions and the y axis the modularity
observed with that resolution.

## Usage

``` r
plot_resolution_res(object, print_head = TRUE, ...)
```

## Arguments

- object:

  The class, either `RbhGraph` or `BulkCoExp`.

- print_head:

  Boolean. Print the Top5 resolution parameters and their meta data.
  Only applicable for `BulkCoExp` objects.

- ...:

  Additional arguments to parse to the functions.

## Value

Plots the result, if the results were found in the class. Otherwise,
throws a warning and returns NULL.

## Examples

``` r
# modularity across the resolutions tested on an RBH graph
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
object <- find_rbh_communities(object, parallel = FALSE, .verbose = FALSE)
plot_resolution_res(object)
```
