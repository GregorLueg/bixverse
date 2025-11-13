# Generate an RBH graph.

This function will generate an RBH graph based on set similarity between
gene modules. You have the option to use an overlap coefficient instead
of Jaccard similarity and to specify a minimum similarity.

## Usage

``` r
generate_rbh_graph(
  object,
  minimum_similarity,
  overlap_coefficient = FALSE,
  spearman = FALSE
)
```

## Arguments

- object:

  The underlying class, see [`rbh_graph()`](rbh_graph.md).

- minimum_similarity:

  The minimum similarity to create an edge.

- overlap_coefficient:

  Boolean. Shall the overlap coefficient be used instead of Jaccard
  similarity. Only relevant if the underlying class is set to set
  similarity.

- spearman:

  Boolean. Shall Spearman correlation be used. Only relevant if the
  underlying class is set to correlation-based similarity.

## Value

The class with added properties.
