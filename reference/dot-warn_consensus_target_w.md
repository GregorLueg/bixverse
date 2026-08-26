# Warn when clustering in cell space is about to get expensive

`consensus_target = "w"` clusters in sample space. On bulk that is
cheap, since there are tens of samples. On single cell the samples are
cells, so it pools a dense `(k * n_runs) x n_cells` matrix and runs an
exhaustive cosine search over rows of that width. Worth a heads-up
before someone waits an hour for it.

## Usage

``` r
.warn_consensus_target_w(nmf_consensus_params, n_samples, k, n_runs)
```

## Arguments

- nmf_consensus_params:

  List. Output of
  [`params_nmf_consensus()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_consensus.md).

- n_samples:

  Integer. Number of cells (or meta cells) in the run.

- k:

  Integer. Components per restart.

- n_runs:

  Integer. Number of restarts.

## Value

`NULL`, invisibly. Called for the warning.
