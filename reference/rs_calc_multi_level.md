# Calculates p-values for pre-processed data

Calculates p-values for pre-processed data

## Usage

``` r
rs_calc_multi_level(stats, es, pathway_size, sample_size, seed, eps, sign)
```

## Arguments

- stats:

  Named numerical vector. Needs to be sorted. The gene level statistics.

- es:

  Numerical vector. The enrichment scores of the pathways of that
  specific size

- pathway_size:

  Integer. The size of the pathways to test.

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

- pvals The pvalues.

- is_cp_ge_half Flag indicating if conditional probability is ≥ 0.5.
  Indicates overesimation of the p-values.
