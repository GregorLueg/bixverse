# Wrapper function for WNN parameters

Wrapper function for WNN parameters

## Usage

``` r
params_sc_wnn(
  k_nn = 20L,
  knn_range = 200L,
  sigma_method = c("snn_farthest", "sigma_idx"),
  sigma_idx = 19L,
  snn_type = c("full_connection", "limited"),
  s_nn = 20L,
  sd_scale = 1,
  kernel_power = 1,
  cross_const = 1e-04,
  sigma_floor = 1e-08,
  knn = list()
)
```

## Arguments

- k_nn:

  Integer. Final number of multimodal neighbours per cell. Defaults to
  `20L`.

- knn_range:

  Integer. Candidate pool size per modality. Each cell's kNN input must
  contain at least this many neighbours. Defaults to `200L`.

- sigma_method:

  String. Bandwidth method. One of `c("snn_farthest", "sigma_idx")`.
  Defaults to `"snn_farthest"`.

- sigma_idx:

  Integer. `"sigma_idx"` only: 0-based kNN index for bandwidth.e.
  `k_nn - 1`). Defaults to `19L`.

- snn_type:

  String. sNN type. The limited version only considers edges that exist
  in the kNN. One of `c("full_connection", "limited")`. Defaults to
  `"full_connection"`.

- s_nn:

  Integer. `"snn_farthest"` only: kNN size used to build the SNN graph.
  Defaults to `20L`.

- sd_scale:

  Numeric. Multiplier on sigma. Defaults to `1.0`.

- kernel_power:

  Numeric. Kernel exponent power. Defaults to `1.0`.

- cross_const:

  Numeric. Cross-modality kernel stabiliser. Defaults to `1e-04`.

- sigma_floor:

  Numeric. Minimum sigma value (avoids division by zero). Defaults to
  `1e-08`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for available parameters: `k`, `knn_method`, `ann_dist`,
  `search_budget`, `n_trees`, `delta`, `diversify_prob`, `ef_budget`,
  `extract_knn`, `m`, `ef_construction`, `ef_search`, `n_list` and
  `n_probe`. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for the available elements. Defaults to
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A named list with the following elements:

- k_nn - Integer. Final number of multimodal neighbours per cell.
  Defaults to `20L`.

- knn_range - Integer. Candidate pool size per modality. Each cell's kNN
  input must contain at least this many neighbours. Defaults to `200L`.

- sigma_method - String. Bandwidth method. One of
  `c("snn_farthest", "sigma_idx")`. Defaults to `"snn_farthest"`.

- sigma_idx - Integer. `"sigma_idx"` only: 0-based kNN index for
  bandwidth.e. `k_nn - 1`). Defaults to `19L`.

- snn_type - String. sNN type. The limited version only considers edges
  that exist in the kNN. One of `c("full_connection", "limited")`.
  Defaults to `"full_connection"`.

- s_nn - Integer. `"snn_farthest"` only: kNN size used to build the SNN
  graph. Defaults to `20L`.

- sd_scale - Numeric. Multiplier on sigma. Defaults to `1.0`.

- kernel_power - Numeric. Kernel exponent power. Defaults to `1.0`.

- cross_const - Numeric. Cross-modality kernel stabiliser. Defaults to
  `1e-04`.

- sigma_floor - Numeric. Minimum sigma value (avoids division by zero).
  Defaults to `1e-08`.

- The elements of
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md),
  overridden by `knn`, spliced in at this position.
