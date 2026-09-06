# Plot the LDA topic count sweep

The three selection metrics and the combined score against `k`. Look for
the `k` where the combined score peaks, then sanity check that Cao Juan
has not already started climbing, which is the signal that topics are
duplicating.

## Usage

``` r
# S3 method for class 'LdaKSweepResult'
plot(x, ...)
```

## Arguments

- x:

  `LdaKSweepResult` object.

- ...:

  Ignored.

## Value

A `ggplot2` object with one panel per metric.

## Examples

``` r
# the three selection metrics and the combined score against k
set.seed(42L)
corpus <- matrix(rbinom(200L * 40L, 1L, 0.05), nrow = 200L, ncol = 40L)
corpus[1:100, 1:10] <- rbinom(1000L, 1L, 0.6)
corpus[101:200, 11:20] <- rbinom(1000L, 1L, 0.6)
colnames(corpus) <- sprintf("term_%02d", 1:40)
sweep_res <- lda_k_sweep(corpus > 0, k_range = 5:7, .verbose = FALSE)
plot(sweep_res)
```
