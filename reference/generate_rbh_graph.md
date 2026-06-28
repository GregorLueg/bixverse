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
