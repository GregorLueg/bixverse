# Print a blitzGSEA null model

Print a blitzGSEA null model

## Usage

``` r
# S3 method for class 'BlitzGseaNull'
print(x, ...)
```

## Arguments

- x:

  `BlitzGseaNull` object.

- ...:

  Ignored.

## Value

Invisibly returns `x`.

## Examples

``` r
# the calibration summary of a blitzGSEA null
set.seed(42L)
stats <- stats::setNames(rnorm(500), sprintf("gene_%03i", 1:500))
print(blitzgsea_calibrate(
  stats,
  blitz_params = params_blitzgsea(permutations = 1000L, anchors = 10L)
))
#> BlitzGseaNull (calibrated blitzGSEA null model)
#>   Signature:        500 genes
#>   Anchors:          10 (sizes 1 to 250)
#>   Centred:          TRUE
#>   KS p-value:       0.404 positive tail, 0.297 negative tail
```
