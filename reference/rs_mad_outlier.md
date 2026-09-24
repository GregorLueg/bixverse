# Calculate MAD outlier detection in Rust.

**\[experimental\]**

## Usage

``` r
rs_mad_outlier(x, threshold, direction)
```

## Arguments

- x:

  Numerical vector to test.

- threshold:

  Numeric. Number of (unscaled) MADs from the median in either direction
  that is acceptable.

- direction:

  String. One of `c("below", "above", "twosided")`. Shall the outlier
  direction be done for values below the threshold, above the threshold
  or in both directions. Unknown strings default to twosided tests.

## Value

A list with the following items:

- outlier - Boolean vector if element is an outlier

- threshold - Applied margin, i.e. `threshold * MAD`.

## Details

Should you provide an empty vector, the function will return an empty
boolean and a threshold of 0.
