# Sweep k and report consensus stability against reconstruction error

**\[experimental\]** Returns diagnostics only, no factors, so a wide
`k_range` stays cheap in memory. The counts are loaded once and reused
across every k. Pick the k where stability is high and the error curve
has not yet flattened, then call
[`rs_nmf_consensus_sc()`](https://gregorlueg.github.io/bixverse/reference/rs_nmf_consensus_sc.md)
there.

## Usage

``` r
rs_nmf_k_sweep_sc(
  f_path_gene,
  gene_indices,
  cell_indices,
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

- f_path_gene:

  Path to the `counts_genes.bin` file.

- gene_indices:

  Integer vector. 0-indexed(!) positions of the genes to include.

- cell_indices:

  Integer vector. 0-indexed(!) positions of cells to include in the
  analysis.

- k_range:

  Integer vector. Ranks to evaluate, every entry at least 2.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- use_second_layer:

  Boolean. If `TRUE`, runs NMF on the normalised counts; if `FALSE`, on
  the raw counts.

- nmf_hals_params:

  Named list. Contains the NMF parameters.

- nmf_consensus_params:

  Named list. Contains the consensus parameters.

- n_runs:

  Integer. Number of restarts per k. Must be at least 2.

- seed:

  Integer. Base random seed.

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
