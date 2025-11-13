# Calculate MAD outlier detection in Rust.

Calculate MAD outlier detection in Rust.

## Usage

``` r
rs_mad_outlier(x, threshold, direction)
```

## Arguments

- x:

  Numerical vector to test.

- threshold:

  Numeric. Number of MADs in either direction that is acceptable.

- direction:

  String. One of `c("below", "above", "twosided")`. Shall the outlier
  direction be done for values below the threshold, above the threshold
  or in both directions.

## Value

A list with the following items:

- outlier - Boolean vector if element is an outlier

- threshold - Applied final threshold

## Details

Should you provide too short vectors, the function will return an empty
boolean and a threshold of 0.
