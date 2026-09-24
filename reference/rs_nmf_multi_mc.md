# Run multiple NMF (HALS) restarts on MetaCells

**\[experimental\]** Assumes that the sparse data is pre-filtered for
the cells/genes you wish to include. Indices in the sparse data need to
be 0-indexed. Both data layers hold the supplied values, so which assay
NMF runs on is decided by what is passed in, not by `use_second_layer`.

## Usage

``` r
rs_nmf_multi_mc(
  sparse_data,
  k,
  preprocessing,
  use_second_layer,
  nmf_hals_params,
  n_runs,
  seed,
  verbose
)
```

## Arguments

- sparse_data:

  A named list with `data`, `indptr`, `indices`, `nrow`, `ncol` and
  `cs_type`. Shape is (metacells, genes).

- k:

  Integer. Number of latent factors per run.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- use_second_layer:

  Boolean. Shall the second data layer be used.

- nmf_hals_params:

  Named list. Contains the NMF parameters, see
  [`params_nmf_hals()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_hals.md).
  The `nmf_init` field is ignored, restarts always use random
  initialisation.

- n_runs:

  Integer. Number of random restarts.

- seed:

  Integer. Base random seed. Run `i` uses `seed + i`.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- w_all - Column-bound `W` matrices across all runs, shape
  `n_meta_cells x (k * n_runs)`.

- h_per_run - List of `H` matrices, each `k x n_genes`.

- losses - Numeric vector. Final reconstruction loss per run.

- converged - Logical vector. Convergence flag per run.

- best_idx - Integer. 1-indexed position of the run with the lowest
  final loss.
