# Generate a new NEBULA result

Wraps the NEBULA fits together. The per-gene test sits in `results`, the
full fixed-effect table in `coefficients` and `se`, so the variance
decomposition is available without re-running anything.

## Usage

``` r
new_sc_nebula_res(results, coefficients, se, params)
```

## Arguments

- results:

  data.table. One row per gene that survived NEBULA's own expression
  filter.

- coefficients:

  Numeric matrix of genes x coefficients. The fixed effects on the
  design scale.

- se:

  Numeric matrix of genes x coefficients. The matching standard errors.

- params:

  Named list. The parameters that were used to generate these results.

## Value

A `ScNebula` class holding the provided data.

## References

He, et al., Commun Biol, 2021
