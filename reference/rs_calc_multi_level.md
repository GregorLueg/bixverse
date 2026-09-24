# Calculates p-values for pre-processed data

**\[experimental\]**

## Usage

``` r
rs_calc_multi_level(stats, es, pathway_size, sample_size, seed, eps, sign)
```

## Arguments

- stats:

  Named numerical vector. Needs to be sorted. The gene level statistics.

- es:

  Numerical vector. The enrichment scores of the pathways.

- pathway_size:

  Integer vector. The size of each pathway, same length as `es`.

- sample_size:

  Integer. The size of the random gene sets to test against.

- seed:

  Integer. Random seed.

- eps:

  Float. Boundary for calculating the p-value.

- sign:

  Boolean. Used for the only positive or only negative score version.

## Value

List with the following elements:

- pvals The p-values.

- is_cp_ge_half Flag indicating if conditional probability is `>= 0.5`.
  Indicates overestimation of the p-values.
