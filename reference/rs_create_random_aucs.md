# Create random AUCs

This function creates a random set of AUCs based on a score vector and a
size of the positive set. This can be used for permutation- based
estimation of Z-scores and subsequently p-values.

## Usage

``` r
rs_create_random_aucs(score_vec, size_pos, random_iters, auc_iters, seed)
```

## Arguments

- score_vec:

  The overall vector of scores.

- size_pos:

  The size of the hits represented in the score_vec.

- random_iters:

  Number of random AUCs to generate.

- auc_iters:

  Number of random iterations to approximate the AUCs. Recommended size:
  10000L.

- seed:

  Seed.

## Value

A vector of random AUCs based the score vector and size of the positive
set.
