# Constructor for MELD parameters

Constructor for MELD parameters

## Usage

``` r
params_meld(
  beta = 60,
  offset = 0,
  order = 1,
  filter = c("heat", "laplacian"),
  chebyshev_order = 50L,
  lap_type = c("combinatorial", "normalised"),
  normalise_indicators = TRUE,
  knn = list()
)
```

## Arguments

- beta:

  Numeric. Smoothing strength; larger values produce smoother densities.
  Must be strictly positive. Defaults to `60.0`.

- offset:

  Numeric. Shift of the filter centre in the rescaled spectrum. Must be
  in `[0, 1]`. Defaults to `0.0`.

- order:

  Numeric. Filter falloff sharpness; larger values approach a square
  low-pass. Must be strictly positive. Defaults to `1.0`.

- filter:

  String. Filter family to use. One of `c("heat", "laplacian")`.
  Defaults to `"heat"`.

- chebyshev_order:

  Integer. Number of Chebyshev coefficients (polynomial terms). Must be
  \>= 2. Defaults to `50L`.

- lap_type:

  String. Type of Laplacian to use for spectral filtering. One of
  `c("combinatorial", "normalised")`. Defaults to `"combinatorial"`.

- normalise_indicators:

  Boolean. If `TRUE`, each column of the indicator matrix is divided by
  its column sum before filtering, making cross-condition densities
  comparable regardless of cells-per-condition. Defaults to `TRUE`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for available parameters: `k`, `knn_method`, `ann_dist`,
  `search_budget`, `n_trees`, `delta`, `diversify_prob`, `ef_budget`,
  `m`, `ef_construction`, `ef_search`, `n_list` and `n_probe`. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for the available elements. Defaults to
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A named list with the following elements:

- The elements of
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md),
  overridden by `knn`, spliced in at this position.

- beta - Numeric. Smoothing strength; larger values produce smoother
  densities. Must be strictly positive. Defaults to `60.0`.

- offset - Numeric. Shift of the filter centre in the rescaled spectrum.
  Must be in `[0, 1]`. Defaults to `0.0`.

- order - Numeric. Filter falloff sharpness; larger values approach a
  square low-pass. Must be strictly positive. Defaults to `1.0`.

- filter - String. Filter family to use. One of
  `c("heat", "laplacian")`. Defaults to `"heat"`.

- chebyshev_order - Integer. Number of Chebyshev coefficients
  (polynomial terms). Must be \>= 2. Defaults to `50L`.

- lap_type - String. Type of Laplacian to use for spectral filtering.
  One of `c("combinatorial", "normalised")`. Defaults to
  `"combinatorial"`.

- normalise_indicators - Boolean. If `TRUE`, each column of the
  indicator matrix is divided by its column sum before filtering, making
  cross-condition densities comparable regardless of
  cells-per-condition. Defaults to `TRUE`.
