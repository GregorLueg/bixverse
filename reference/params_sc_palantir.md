# Wrapper function for Palantir parameters

Parameters controlling the Palantir trajectory inference. The kNN graph
you hand to
[`run_palantir_sc()`](https://gregorlueg.github.io/bixverse/reference/run_palantir_sc.md)
feeds the diffusion kernel. The geodesics are measured over a second kNN
graph that Palantir builds internally on the multiscale space, and `knn`
is what controls that one. It is a different knob from `k` in the kNN
parameter block, which sizes the backend index. Palantir overrides `k`
and `ann_dist` for its internal search, so only `knn_method` and the
backend tuning parameters have an effect.

## Usage

``` r
params_sc_palantir(
  n_dcs = 10L,
  n_eigs = NULL,
  knn = 30L,
  num_waypoints = 1200L,
  scale_components = TRUE,
  use_early_cell_as_start = TRUE,
  max_iterations = 25L,
  branch_prob_threshold = 0.01,
  lanczos_basis_size = NULL,
  lanczos_max_restarts = 16L,
  lanczos_tol = 1e-08,
  knn_params = list()
)
```

## Arguments

- n_dcs:

  Integer. Diffusion components to extract before the multiscale
  scaling. Defaults to `10L`.

- n_eigs:

  Integer or `NULL`. Eigenvectors to retain, not components: the trivial
  leading eigenvector is counted here and then dropped, so `3L` leaves
  two multiscale components. If `NULL`, the count is picked from the
  largest eigengap, as the reference does. Defaults to `NULL`.

- knn:

  Integer. Neighbours for the geodesic graph over the multiscale space,
  in the reference's self-inclusive convention. Defaults to `30L`.

- num_waypoints:

  Integer. Target waypoint count for the max-min sampler. Defaults to
  `1200L`.

- scale_components:

  Boolean. Min-max scale each multiscale component to `[0, 1]` before
  any distance is taken. Defaults to `TRUE`.

- use_early_cell_as_start:

  Boolean. Use the provided early cell directly rather than snapping it
  to the nearest diffusion-map boundary cell. Defaults to `TRUE`.

- max_iterations:

  Integer. Iteration cap for the pseudotime refinement. Defaults to
  `25L`.

- branch_prob_threshold:

  Numeric. Fate probabilities below this are zeroed. Defaults to `0.01`.

- lanczos_basis_size:

  Integer or `NULL`. Krylov basis vectors held at once during the
  diffusion-map eigendecomposition. If `NULL`, derived from the
  requested component count. Defaults to `NULL`.

- lanczos_max_restarts:

  Integer. Maximum restart cycles for the Lanczos solver. Defaults to
  `16L`.

- lanczos_tol:

  Numeric. Relative residual tolerance for the Lanczos solver. Defaults
  to `1e-08`.

- knn_params:

  List. Optional overrides for the kNN parameters of the internal
  multiscale search. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for available parameters: `k`, `knn_method`, `ann_dist`,
  `search_budget`, `n_trees`, `delta`, `diversify_prob`, `ef_budget`,
  `m`, `ef_construction`, `ef_search`, `n_list` and `n_probe`. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for the available elements. Defaults to
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A named list with the following elements:

- n_dcs - Integer. Diffusion components to extract before the multiscale
  scaling. Defaults to `10L`.

- n_eigs - Integer or `NULL`. Eigenvectors to retain, not components:
  the trivial leading eigenvector is counted here and then dropped, so
  `3L` leaves two multiscale components. If `NULL`, the count is picked
  from the largest eigengap, as the reference does. Defaults to `NULL`.

- knn - Integer. Neighbours for the geodesic graph over the multiscale
  space, in the reference's self-inclusive convention. Defaults to
  `30L`.

- num_waypoints - Integer. Target waypoint count for the max-min
  sampler. Defaults to `1200L`.

- scale_components - Boolean. Min-max scale each multiscale component to
  `[0, 1]` before any distance is taken. Defaults to `TRUE`.

- use_early_cell_as_start - Boolean. Use the provided early cell
  directly rather than snapping it to the nearest diffusion-map boundary
  cell. Defaults to `TRUE`.

- max_iterations - Integer. Iteration cap for the pseudotime refinement.
  Defaults to `25L`.

- branch_prob_threshold - Numeric. Fate probabilities below this are
  zeroed. Defaults to `0.01`.

- lanczos_basis_size - Integer or `NULL`. Krylov basis vectors held at
  once during the diffusion-map eigendecomposition. If `NULL`, derived
  from the requested component count. Defaults to `NULL`.

- lanczos_max_restarts - Integer. Maximum restart cycles for the Lanczos
  solver. Defaults to `16L`.

- lanczos_tol - Numeric. Relative residual tolerance for the Lanczos
  solver. Defaults to `1e-08`.

- The elements of
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md),
  overridden by `knn_params`, spliced in at this position.

## References

Setty, et al., Nat. Biotechnol., 2019.
