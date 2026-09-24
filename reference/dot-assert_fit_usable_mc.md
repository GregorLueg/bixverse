# Fetch a meta cell residual fit that is safe to compute with

Meta cells have no cell filter to drift, their cell set is fixed at
construction, so this is the presence check alone.

## Usage

``` r
.assert_fit_usable_mc(object)
```

## Arguments

- object:

  `MetaCells` class.

## Value

The `ScResidualFit`.
