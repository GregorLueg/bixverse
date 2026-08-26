# Run consensus NMF over a set of single cells and genes

**\[experimental\]** Runs `n_runs` HALS NMF restarts, pools their
components, drops unstable ones by local density, k-means clusters the
survivors and refits the partner factor against the per-cluster median.
The counts are streamed from the binary file, but the restart factors
are dense and all held at once, so `n_runs` times `k` times the cell
count is the memory you should budget for.

## Usage

``` r
rs_nmf_consensus_sc(
  f_path_gene,
  gene_indices,
  cell_indices,
  k,
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

- k:

  Integer. Number of latent factors. Must be at least 2.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- use_second_layer:

  Boolean. If `TRUE`, runs NMF on the normalised counts; if `FALSE`, on
  the raw counts.

- nmf_hals_params:

  Named list. Contains the NMF parameters. The `nmf_init` field is
  ignored, restarts always use random initialisation.

- nmf_consensus_params:

  Named list. Contains the consensus parameters.

- n_runs:

  Integer. Number of restarts. Must be at least 2.

- seed:

  Integer. Base random seed. Restart `i` uses `seed + i`.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- w - The left factor matrix (n_cells x k)

- h - The right factor matrix (k x n_genes)

- rel_error - Reconstruction error relative to the squared Frobenius
  norm of the input. Not comparable with the absolute `final_loss` the
  single-run version returns.

- rel_run_errors - The same, per restart.

- labels - Integer vector of length `k * n_runs`. Cluster each pooled
  component landed in, `NA` if it was dropped.

- local_density - Mean cosine distance to the nearest neighbours per
  pooled component.

- kept - 1-indexed positions of the surviving pooled components.

- silhouette - Silhouette per survivor, aligned with `kept`.

- stability - Mean silhouette over the survivors.

- cluster_sizes - Number of survivors per cluster.

- n_dropped - Number of pooled components removed.

- n_empty_clusters - Number of clusters left with no members.

## References

Kotliar et al., eLife, 2019
