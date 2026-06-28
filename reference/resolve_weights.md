# Resolve prioritisation weights from the scenario or validate a user vector

Resolve prioritisation weights from the scenario or validate a user
vector

## Usage

``` r
resolve_weights(weights = NULL, scenario, has_condition_de)
```

## Arguments

- weights:

  Optional named numeric.

- scenario:

  String. One of `c("case_control", "one_condition")`.

- has_condition_de:

  Boolean. If `scenario = "case_control"`, validates that the
  corresponding weights are non-zero.

## Value

Final weights
