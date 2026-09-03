# Run Palantir trajectory inference over a kNN graph

**\[experimental\]** Implementation of Palantir in Rust. The provided
kNN graph feeds the diffusion kernel. The geodesics themselves are
measured over a second kNN graph that is built internally on the
multiscale space, which is where the reference measures them.

## Usage

``` r
rs_palantir(
  knn_data,
  palantir_params,
  early_cell,
  terminal_states,
  seed,
  verbose
)
```

## Arguments

- knn_data:

  List. The `SingleCellNearestNeighbour` data with `indices`
  (0-indexed!), `dist`, `k` and `dist_metric`.

- palantir_params:

  List. Parameter list, see
  [`params_sc_palantir()`](https://gregorlueg.github.io/bixverse/reference/params_sc_palantir.md).

- early_cell:

  Integer. Index (0-indexed!) of the early cell within the rows of the
  kNN data.

- terminal_states:

  Optional integer vector. Indices (0-indexed!) of the terminal states.
  If `NULL`, they are detected from the waypoint Markov chain.

- seed:

  Integer. Seed for reproducibility purposes.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with

- pseudotime - Numerical vector with the pseudotime per cell, min-max
  scaled to `[0, 1]`. The start cell is not pinned to 0; a start cell
  far from 0 means the refinement disagreed with the anchor.

- entropy - Numerical vector with the differentiation entropy per cell
  (natural log).

- branch_probs - Numerical matrix of cells x terminal states with the
  fate probabilities. Rows need not sum to one, as sub-threshold values
  are zeroed without renormalisation.

- terminal_states - Integer vector with the terminal state cell indices
  (0-indexed!). Sets the column order of `branch_probs`.

- waypoints - Integer vector with the waypoint cell indices
  (0-indexed!). The first element is the start cell.

- start_cell - Integer. The start cell that was actually used
  (0-indexed!).

- multiscale - Numerical matrix of cells x components with the
  multiscale diffusion components.

- iterations - Integer. Refinement passes that were run.

- converged - Boolean. Did the pseudotime refinement converge before the
  cap.

- eigen_converged - Boolean. Did the diffusion eigensolve meet its
  tolerance rather than running out of restarts. `FALSE` means the
  embedding is under-resolved and every distance taken on it is suspect.

- eigen_residual - Numeric. Largest achieved `||A x - lambda x||` from
  the diffusion eigensolve.

- repair_edges - Integer. Bridging edges the connectivity repair had to
  add. Anything non-zero means the kNN graph was disconnected.

- stranded_waypoints - Integer. Waypoints from which no terminal state
  is reachable.

## References

Setty, et al., Nat. Biotechnol., 2019.
