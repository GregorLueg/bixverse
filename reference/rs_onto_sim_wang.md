# Calculate the Wang similarity for specific terms

This function calculates the Wang similarities between all permutations
of a given set of terms.

## Usage

``` r
rs_onto_sim_wang(terms, parents, children, w)
```

## Arguments

- terms:

  String vector. The terms you wish to calculate the similarities for.

- parents:

  String vector. The names of the parents.

- children:

  String vector. The names of the childs. The length of `parents` needs
  to be equal to `children`.

- w:

  Numerics. The weights between the parents and children. Need to be
  values between 0 and 1.

## Value

A list with:

- term1 - name of the first term.

- term2 - name of the second term.

- sims - similarity between the two terms.
