# Generates a new `ADTCounts` class via DSB normalisation

This function generates a new `ADTCounts` class using DSB normalisation
from Mulè et al., instead of CLR. When `empty_drops` is provided, the
per-protein background is estimated from empty droplets. Without
`empty_drops`, the function falls back to a 2-component k-means on the
log-transformed cell counts ("ModelNegativeADTnorm").

## Usage

``` r
new_adt_counts_dsb(
  raw_counts,
  cell_info,
  empty_drops = NULL,
  isotype_names = NULL,
  dsb_params = params_sc_dsb(),
  scale_factor = c("standardise", "mean_subtract"),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- raw_counts:

  Numeric matrix. Cells x proteins matrix of raw ADT counts with cell
  barcodes as row names and protein names as column names.

- cell_info:

  Named integer vector. Output of
  [`get_cell_info()`](https://gregorlueg.github.io/bixverse/reference/get_cell_info.md).
  Defines as elements the cell indices (R-based) and as names the
  barcodes.

- empty_drops:

  Optional numeric matrix. Cells x proteins matrix of empty-droplet ADT
  counts. If provided, used to estimate per-protein ambient background.

- isotype_names:

  Optional string vector. Column names in `raw_counts` identifying
  isotype control proteins. Required when
  `dsb_params$use_isotype_controls = TRUE`.

- dsb_params:

  List. Output of
  [`params_sc_dsb()`](https://gregorlueg.github.io/bixverse/reference/params_sc_dsb.md)
  with DSB parameters.

- scale_factor:

  String. One of `c("standardise", "mean_subtract")`. Only used when
  `empty_drops` is provided.

- seed:

  Integer. Random seed for k-means initialisation.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

`ADTCounts` containing the raw counts and the DSB-normalised counts.

## References

Mulè et al., Nat Commun, 2022
