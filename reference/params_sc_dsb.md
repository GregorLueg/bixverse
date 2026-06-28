# Default parameters for DSB ADT normalisation

Default parameters for DSB ADT normalisation

## Usage

``` r
params_sc_dsb(
  denoise_counts = TRUE,
  use_isotype_controls = TRUE,
  pseudocount = 10,
  quantile_low = NULL,
  quantile_high = NULL
)
```

## Arguments

- denoise_counts:

  Boolean. Run Step II (cell-to-cell technical noise removal).

- use_isotype_controls:

  Boolean. Include isotype controls in the noise matrix in Step II.
  Requires `isotype_indices` to be passed at call time.

- pseudocount:

  Numeric. Pseudocount added before the log transform. The DSB paper
  recommends `10` with empty droplets and `1` without.

- quantile_low:

  Optional numeric in `[0, 1)`. Lower quantile for per-protein output
  clipping. If `NULL` (and `quantile_high` is also `NULL`), no clipping
  is applied.

- quantile_high:

  Optional numeric in `(0, 1]`. Upper quantile for per-protein output
  clipping. If `NULL` (and `quantile_low` is also `NULL`), no clipping
  is applied.

## Value

A list with the parameters.
