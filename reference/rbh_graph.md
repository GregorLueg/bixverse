# Reciprocal best hit graph (deprecated)

**\[deprecated\]**

This constructor has been renamed to
[`RbhGraph()`](https://gregorlueg.github.io/bixverse/reference/RbhGraph.md).

## Usage

``` r
rbh_graph(
  module_results,
  rbh_type = c("cor", "set"),
  dataset_col = NULL,
  module_col = NULL,
  value_col = NULL
)
```

## Arguments

- module_results:

  data.table or list containing the gene modules.

- rbh_type:

  String. One of `c("cor", "set")`.

- dataset_col:

  The column indicating the data set/method origin.

- module_col:

  The column storing the module names.

- value_col:

  The column storing the genes within each module.

## Value

Returns a
[`RbhGraph()`](https://gregorlueg.github.io/bixverse/reference/RbhGraph.md)
object.
