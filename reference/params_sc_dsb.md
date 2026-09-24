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

  Boolean. Run Step II (cell-to-cell technical noise removal). Defaults
  to `TRUE`.

- use_isotype_controls:

  Boolean. Include isotype controls in the noise matrix in Step II.
  Requires `isotype_indices` to be passed at call time. Defaults to
  `TRUE`.

- pseudocount:

  Numeric. Pseudocount added before the log transform. The DSB paper
  recommends `10` with empty droplets and `1` without. Defaults to
  `10.0`.

- quantile_low:

  Numeric or `NULL`. Optional numeric in `[0, 1)`. Lower quantile for
  per-protein output clipping. If `NULL` (and `quantile_high` is also
  `NULL`), no clipping is applied. Defaults to `NULL`.

- quantile_high:

  Numeric or `NULL`. Optional numeric in `(0, 1]`. Upper quantile for
  per-protein output clipping. If `NULL` (and `quantile_low` is also
  `NULL`), no clipping is applied. Defaults to `NULL`.

## Value

A named list with the following elements:

- denoise_counts - Boolean. Run Step II (cell-to-cell technical noise
  removal). Defaults to `TRUE`.

- use_isotype_controls - Boolean. Include isotype controls in the noise
  matrix in Step II. Requires `isotype_indices` to be passed at call
  time. Defaults to `TRUE`.

- pseudocount - Numeric. Pseudocount added before the log transform. The
  DSB paper recommends `10` with empty droplets and `1` without.
  Defaults to `10.0`.

- quantile_low - Numeric or `NULL`. Optional numeric in `[0, 1)`. Lower
  quantile for per-protein output clipping. If `NULL` (and
  `quantile_high` is also `NULL`), no clipping is applied. Defaults to
  `NULL`.

- quantile_high - Numeric or `NULL`. Optional numeric in `(0, 1]`. Upper
  quantile for per-protein output clipping. If `NULL` (and
  `quantile_low` is also `NULL`), no clipping is applied. Defaults to
  `NULL`.

\[0,
1)`. Lower quantile for per-protein output clipping. If `NULL`(and`quantile_high`is also`NULL`), no clipping is applied. Defaults to `NULL`. \item quantile_high - Numeric or `NULL`. Optional numeric in `(0,
1\]:
R:0,%201)%60.%20Lower%0A%20quantile%20for%20per-protein%20output%20clipping.%20If%20%60NULL%60%20(and%20%60quantile_high%60%20is%0A%20also%20%60NULL%60),%20no%20clipping%20is%20applied.%20Defaults%20to%20%60NULL%60.%0A%20%5C%5Citem%20quantile_high%20-%20Numeric%20or%20%60NULL%60.%20Optional%20numeric%20in%20%60(0,%201
