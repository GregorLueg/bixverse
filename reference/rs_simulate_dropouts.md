# Sparsify bulkRNAseq like data

This function takes in a (raw) count matrix (for example from the
synthetic data in bixverse) and applies sparsification to it based on
two possible functions:

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
rs_simulate_dropouts(
  count_mat,
  dropout_function,
  dropout_midpoint,
  dropout_shape,
  power_factor,
  global_sparsity,
  seed
)
```

## Arguments

- count_mat:

  Numerical matrix. Original numeric matrix.

- dropout_function:

  String. One of `c("log", "powerdecay")`. Defines which function will
  be used to induce the sparsity.

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

  Integer. Seed for reproducibility.

## Value

The sparsified matrix based on the provided parameters.
