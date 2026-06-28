# Default parameters for Symphony query mapping

Default parameters for Symphony query mapping

## Usage

``` r
params_symphony_map(sigma = 0.1, lambda = 1)
```

## Arguments

- sigma:

  Numeric. Soft-clustering fuzziness for query -\> reference centroid
  assignment. Symphony R default is 0.1.

- lambda:

  Numeric. Ridge penalty on batch coefficients. Symphony R hardcodes
  1.0.

## Value

A list with the parameters.
