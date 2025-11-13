# Simulate dropouts via different functions on synthetic bulk data

This function induces expression-level dependent sparsity on the data.
The two possible functions are:

**Logistic function:**

With dropout probability defined as:

    P(dropout) = clamp(1 / (1 + exp(shape * (ln(exp+1) - ln(midpoint+1)))), 0.3, 0.8) * (1 - global_sparsity) + global_sparsity

with the following characteristics:

- Plateaus at global_sparsity dropout for high expression genes

- Partial dropout preserves count structure via binomial thinning

- Good for preserving variance-mean relationships

**Power Decay function:**

With dropout probability defined as:

    P(dropout) = (midpoint / (exp + midpoint))^power * scale_factor * (1 - global_sparsity) + global_sparsity

with the following characteristics:

- No plateau - high expression genes get substantial dropout

- Complete dropout only (no partial dropout)

- More uniform dropout across expression range

## Usage

``` r
simulate_dropouts(
  object,
  dropout_function = c("log", "powerdecay"),
  dropout_midpoint = 5,
  dropout_shape = 0.5,
  power_factor = 0.3,
  global_sparsity = 0.25,
  seed = 123L
)
```

## Arguments

- object:

  The `synthetic_bulk_data` class.

- dropout_function:

  String. One of `c("log", "powerdecay")`. These are the two options to
  cause dropout in the bulk data.

- dropout_midpoint:

  Numeric. Controls the midpoint parameter of the logistic and power
  decay function.

- dropout_shape:

  Numeric. Controls the shape parameter of the logistic function.

- power_factor:

  Numeric. Controls the power factor of the power decay function.

- global_sparsity:

  Numeric. The global sparsity parameter.

- seed:

  Integer. Random seed for reproducibility purposes.

## Value

`synthetic_bulk_data` with added sparse data.
