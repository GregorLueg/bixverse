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

## Examples

``` r
# the five highest probability terms of each topic
set.seed(42L)
corpus <- matrix(rbinom(200L * 40L, 1L, 0.05), nrow = 200L, ncol = 40L)
corpus[1:100, 1:10] <- rbinom(1000L, 1L, 0.6)
corpus[101:200, 11:20] <- rbinom(1000L, 1L, 0.6)
colnames(corpus) <- sprintf("term_%02d", 1:40)
lda_res <- run_lda(corpus > 0, k = 2L, .verbose = FALSE)
get_top_terms(lda_res, n = 5L)
#>        topic  rank    term probability
#>       <char> <int>  <char>       <num>
#>  1: topic_01     1 term_07  0.09487961
#>  2: topic_01     2 term_09  0.09487961
#>  3: topic_01     3 term_03  0.09221070
#>  4: topic_01     4 term_02  0.09087624
#>  5: topic_01     5 term_10  0.08687288
#>  6: topic_02     1 term_18  0.10072022
#>  7: topic_02     2 term_14  0.09937908
#>  8: topic_02     3 term_12  0.09267333
#>  9: topic_02     4 term_19  0.08864989
#> 10: topic_02     5 term_17  0.08730874
```
