# Apply a range normalisation on a vector.

Applies a range normalisation on an R vector.

## Usage

``` r
rs_range_norm(x, max_val, min_val)
```

## Arguments

- x:

  Numerical vector. The data to normalise.

- max_val:

  Numeric. The upper bound value to normalise into. If set to 1, the
  function will be equal to a min-max normalisation.

- min_val:

  Numeric. The lower bound value to normalise into. If set to 0, the
  function will equal a min-max normalisation.

## Value

Normalised values
