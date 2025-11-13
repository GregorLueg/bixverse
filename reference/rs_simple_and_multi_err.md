# Calculates the simple and multi error for fgsea multi level

Calculates the simple and multi error for fgsea multi level

## Usage

``` r
rs_simple_and_multi_err(n_more_extreme, nperm, sample_size)
```

## Arguments

- n_more_extreme:

  Integer vector. The number of times the ES was larger than the
  permutations.

- nperm:

  Integer. Number of permutations.

- sample_size:

  Integer. Number of samples.

## Value

List with the following elements:

- simple_err Vector of simple errors.

- multi_err Vector of multi errors.
