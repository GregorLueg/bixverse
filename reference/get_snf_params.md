# Get the SNF params

Get the SNF params

## Usage

``` r
get_snf_params(object)
```

## Arguments

- object:

  The underlying class
  [`SimilarityNetworkFusion()`](https://gregorlueg.github.io/bixverse/reference/SimilarityNetworkFusion.md).

## Value

Returns the stored SNF params

## Examples

``` r
# the SNF parameters stored in an empty class
object <- SimilarityNetworkFusion(snf_params = params_snf(k = 3L))
get_snf_params(object)$k
#> [1] 3
```
