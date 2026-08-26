# Get the highest-probability terms per topic

The interpretability surface of a topic model: for each topic, the terms
carrying the most probability mass. Note that these are probabilities
within a topic, so they are comparable down a topic but not across
topics of different breadth.

## Usage

``` r
get_top_terms(x, n = 20L)

# S3 method for class 'LdaResult'
get_top_terms(x, n = 20L)
```

## Arguments

- x:

  `LdaResult` object.

- n:

  Integer. Number of terms to return per topic.

## Value

A data.table with `topic`, `rank`, `term` and `probability`, sorted by
topic and then rank.
