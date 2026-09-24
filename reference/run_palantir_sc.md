# Run Palantir trajectory inference

Palantir models differentiation as a Markov chain over a diffusion map
of the cells. It returns a pseudotime per cell, the probability of each
cell reaching each terminal state, and the differentiation entropy
derived from those probabilities. For details, please refer to Setty, et
al.

The kNN graph stored on the object feeds the diffusion kernel, so run
[`find_neighbours_sc()`](https://gregorlueg.github.io/bixverse/reference/find_neighbours_sc.md)
first. The geodesics themselves are measured over a second kNN graph
that Palantir builds internally on the multiscale space, controlled by
the `knn` element of
[`params_sc_palantir()`](https://gregorlueg.github.io/bixverse/reference/params_sc_palantir.md).

For `SingleCellsMultiModal` the `"wnn"` graph works, but its distances
are kernel-derived and not a metric, so the diffusion map on top is a
heuristic on a heuristic. You will get a warning saying as much.

## Usage

``` r
run_palantir_sc(
  object,
  early_cell,
  terminal_states = NULL,
  modality = c("rna", "adt", "wnn"),
  palantir_params = params_sc_palantir(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  One of `SingleCells`, `SingleCellsSubset`, `MetaCells` or
  `SingleCellsMultiModal`.

- early_cell:

  String. Name of the cell to start the trajectory from. Must be one of
  the cells the kNN graph was built over.

- terminal_states:

  Optional character vector. Names of the terminal state cells. If
  `NULL`, they are detected from the waypoint Markov chain. Name the
  vector (e.g. `c(Ery = "Run4_2005...")`) and those labels become the
  column names of `branch_probs`, which saves relabelling the fate
  columns by hand downstream.

- modality:

  String. One of `c("rna", "adt", "wnn")`. Which kNN graph to run over.
  Anything but `"rna"` requires a `SingleCellsMultiModal` object.

- palantir_params:

  List. See
  [`params_sc_palantir()`](https://gregorlueg.github.io/bixverse/reference/params_sc_palantir.md).

- seed:

  Integer. For reproducibility.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

A `PalantirRes` S3 object with:

- pseudotime - data.table with `cell_id`, `pseudotime` (min-max scaled
  to `[0, 1]`; the start cell is not pinned to 0, and a start cell far
  from 0 means the refinement disagreed with the anchor) and `entropy`.

- branch_probs - Numeric matrix of cells x terminal states with the fate
  probabilities. Rows need not sum to one, as sub-threshold values are
  zeroed without renormalisation.

- terminal_states - Character vector with the terminal state cell names,
  carrying the labels of a named `terminal_states` argument as its own
  names. Sets the column order of `branch_probs`.

- waypoints - Character vector with the waypoint cell names. The first
  element is the start cell.

- start_cell - String. The start cell that was actually used.

- multiscale - Numeric matrix of cells x components with the multiscale
  diffusion components.

- run_info - List with `iterations`, `converged`, `eigen_converged`,
  `eigen_residual`, `repair_edges`, `stranded_waypoints`, `n_waypoints`
  and `modality`. A `eigen_converged` of `FALSE` means the diffusion
  eigensolve ran out of restarts and the embedding is under-resolved.

## References

Setty, et al., Nat. Biotechnol., 2019.

## Examples

``` r
# pseudotime from an arbitrary start cell on the demo kNN graph
sc <- demo_single_cells()
res <- run_palantir_sc(
  sc,
  early_cell = get_knn_obj(sc)$used_cells[1],
  palantir_params = params_sc_palantir(
    knn = 15L,
    num_waypoints = 100L,
    n_eigs = 3L,
    use_early_cell_as_start = TRUE,
    knn_params = list(knn_method = "exhaustive")
  ),
  .verbose = FALSE
)
head(res$pseudotime)
#>     cell_id  pseudotime      entropy
#>      <char>       <num>        <num>
#> 1: cell_001 0.003929917 6.039816e-01
#> 2: cell_002 0.943664372 2.954053e-09
#> 3: cell_003 0.935922384 1.399716e-02
#> 4: cell_004 0.057008285 6.029100e-01
#> 5: cell_005 0.764215469 1.669869e-03
#> 6: cell_006 0.965360820 7.510467e-03

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
