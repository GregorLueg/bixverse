# Wrapper function for HVG detection parameters.

Wrapper function for HVG detection parameters.

## Usage

``` r
params_sc_hvg(
  method = "vst",
  loess_span = 0.3,
  num_bin = 20L,
  bin_method = "equal_width"
)
```

## Arguments

- method:

  String. One of `c("vst", "meanvarbin", "dispersion")`.

- loess_span:

  Numeric. The span parameter for the loess function that is used to
  standardise the variance for `method = "vst"`.

- num_bin:

  Integer. Not yet implemented.

- bin_method:

  String. One of `c("equal_width", "equal_freq")`.

## Value

A list with the HVG parameters
