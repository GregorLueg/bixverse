# Calculates the point wise mutual information

Calculates the pointwise mutual information (can be also normalised)
across all columns of the data. This can be used to identify
(dis)similar samples based on Boolean characteristics.

## Usage

``` r
rs_pointwise_mutual_info(x, normalise)
```

## Arguments

- x:

  Logical matrix. The columns represent features and the rows represent
  samples

- normalise:

  Shall the normalised pointwise mutual information be returned.

## Value

The (normalised) pointwise mutual information matrix.
