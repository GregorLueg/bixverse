# Parameters for the per-cell ScType assignment

Controls the per-cell path of
[`assign_sc_type()`](https://gregorlueg.github.io/bixverse/reference/assign_sc_type.md):
how the raw ScType scores are rescaled, how hard the scores get smoothed
over the sNN graph, and where the cut-offs for an Unknown call and for a
mixed cluster sit.

## Usage

``` r
params_sctype_cells(
  alpha = 0.5,
  iterations = 2L,
  tolerance = 1e-04,
  calibration = c("none", "column_z"),
  score_floor = 0.25,
  purity_threshold = 0.9
)
```

## Arguments

- alpha:

  Numeric in `[0, 1]`. Self-retention during smoothing. Each iteration
  computes `alpha * original + (1 - alpha) * neighbour_average`.

- iterations:

  Integer \>= 0. Number of smoothing iterations. `0` disables smoothing.

- tolerance:

  Numeric \> 0. Convergence tolerance for the smoothing.

- calibration:

  String. One of `c("none", "column_z")`. `"column_z"` standardises each
  cell type's score column across cells, which removes the bias towards
  cell types whose marker sets happen to produce larger scores.

- score_floor:

  Numeric \>= 0. Minimum score for a cell to get a call instead of `NA`.

- purity_threshold:

  Numeric in `[0, 1]`. Cluster purity above which the hybrid assignment
  keeps the cluster-level call.

## Value

A list with the parameters for usage in subsequent functions.

## References

Zhou et al., NIPS, 2004.
