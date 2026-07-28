# Validate hard bounds against MAD direction

Errors if the supplied bounds contradict the metric's MAD direction
(e.g. a `lower` bound on a metric only checked from above).

## Usage

``` r
.check_bounds_direction(bounds, direction, metric_name)
```

## Arguments

- bounds:

  Named numeric vector with `lower` and/or `upper`.

- direction:

  String. One of `"twosided"`, `"below"`, `"above"`.

- metric_name:

  String. Used in the error message.

## Value

Invisible `NULL`. Called for its side effect.
