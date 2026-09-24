# Sweep k and report consensus stability against reconstruction error

**\[experimental\]** Returns diagnostics only, no factors, so a wide
`k_range` stays cheap in memory. Pick the k where stability is high and
the error curve has not yet flattened, then call
[`rs_nmf_consensus_mc()`](https://gregorlueg.github.io/bixverse/reference/rs_nmf_consensus_mc.md)
there. Assumes that the sparse data is pre-filtered for the cells/genes
you wish to include. Indices in the sparse data need to be 0-indexed.
Both data layers hold the supplied values, so which assay NMF runs on is
decided by what is passed in, not by `use_second_layer`.

## Usage

``` r
rs_nmf_k_sweep_mc(
  sparse_data,
  k_range,
  preprocessing,
  use_second_layer,
  nmf_hals_params,
  nmf_consensus_params,
  n_runs,
  seed,
  verbose
)
```

## Arguments

- sparse_data:

  A named list with `data`, `indptr`, `indices`, `nrow`, `ncol` and
  `cs_type`. Shape is (metacells, genes).

- k_range:

  Integer vector. Ranks to evaluate, every entry at least 2.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- use_second_layer:

  Boolean. Shall the second data layer be used.

- nmf_hals_params:

  Named list. Contains the NMF parameters, see
  [`params_nmf_hals()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_hals.md).

- nmf_consensus_params:

  Named list. Contains the consensus parameters, see
  [`params_nmf_consensus()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_consensus.md).

- n_runs:

  Integer. Number of restarts per k. Must be at least 2.

- seed:

  Integer. Base random seed. The i-th k uses `seed + i * n_runs`.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list of equal-length vectors, one element per swept k

- k - The rank.

- stability - Mean silhouette of the consensus clusters. `NaN` where the
  consensus step failed.

- best_error - Lowest restart error, relative to the squared Frobenius
  norm of the input.

- median_error - Median restart error, same scale.

- consensus_failed - Did the density filter leave fewer than `k`
  components.

- n_dropped - Number of pooled components removed.

- n_empty_clusters - Number of clusters left with no members.

- n_converged - Restarts that met the HALS tolerance.

## References

Kotliar et al., eLife, 2019
