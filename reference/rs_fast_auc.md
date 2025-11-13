# Fast AUC calculation

This function calculates rapidly AUCs based on an approximation.

## Usage

``` r
rs_fast_auc(pos_scores, neg_scores, iters, seed)
```

## Arguments

- pos_scores:

  The scores of your hits.

- neg_scores:

  The scores of your non-hits.

- iters:

  Number of iterations to run the function for. Recommended size:
  10000L.

- seed:

  Seed.

## Value

The AUC.
