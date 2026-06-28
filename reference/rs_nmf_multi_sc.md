# Run multiple NMF (HALS) restarts over a set of single cells and genes

**\[experimental\]** Runs `n_runs` HALS NMF with random initialisations
seeded by `seed + i`. The `nmf_init` field in `nmf_hals_params` is
ignored; random init is always used. The returned `w_all` is the
column-binding of all run W matrices.

## Usage

``` r
rs_nmf_multi_sc(
  f_path_gene,
  gene_indices,
  cell_indices,
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

- f_path_gene:

  Path to the `counts_genes.bin` file.

- gene_indices:

  Integer vector. 0-indexed(!) positions of the genes to include.

- cell_indices:

  Integer vector. 0-indexed(!) positions of cells to include in the
  analysis.

- k:

  Integer. Number of latent factors per run.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- use_second_layer:

  Boolean. If `TRUE`, runs NMF on the normalised counts; if `FALSE`, on
  the raw counts.

- nmf_hals_params:

  Named list. Contains the NMF parameters.

- n_runs:

  Integer. Number of random restarts.

- seed:

  Integer. Base random seed. Run `i` uses `seed + i`.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- w_all - Column-bound W matrices across all runs, shape
  `n_features x (k * n_runs)`. Columns `i*k+1..(i+1)*k` are run `i`'s
  components (1-indexed).

- h_per_run - List of H matrices, each `k x n_cells`.

- losses - Numeric vector. Final reconstruction loss per run.

- converged - Logical vector. Convergence flag per run.

- best_idx - Integer. 1-indexed position of the run with the lowest
  final loss.
