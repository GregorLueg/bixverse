# BixverseBaseClass

Generic base class that is used for inheritance in certain common
methods across other classes.

## Usage

``` r
BixverseBaseClass()
```

## Value

Returns the S7 object for further operations.

## Properties

- params:

  A (nested) list that will store all the parameters of the applied
  function.

- final_results:

  A data.table that will contain the final results.

## Examples

``` r
# every analysis class inherits the base class getters
object <- SimilarityNetworkFusion(snf_params = params_snf(k = 3L))
inherits(object, "bixverse::BixverseBaseClass")
#> [1] TRUE
```
