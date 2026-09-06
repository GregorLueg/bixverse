# Calculates a harmonic sum normalised between 0 to 1.

The function takes in a vector of scores between 0 and 1 and calculates
a harmonic sum, based on the approach OpenTargets takes to do their
gene - disease evidence scores, see:
https://platform-docs.opentargets.org/associations.

## Usage

``` r
ot_harmonic_score(x)
```

## Arguments

- x:

  Numeric vector. Needs to be between 0 and 1.

## Value

Harmonic, normalised sum of the provided scores.

## Examples

``` r
# harmonic sum over four evidence scores
ot_harmonic_score(c(1, 0.8, 0.5, 0.1))
#> [1] 0.8863415
```
