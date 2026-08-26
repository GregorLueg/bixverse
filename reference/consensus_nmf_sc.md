# Run consensus NMF on single cell or meta cell data

Runs `n_runs` HALS NMF restarts, pools their components, drops unstable
ones by local density, k-means clusters the survivors into `k` groups
and refits the partner factor against the per-cluster median. This is
cNMF: the answer is the structure the restarts agree on, and the mean
silhouette of those clusters (`stability`) says how much they agreed.

## Usage

``` r
consensus_nmf_sc(
  object,
  k,
  cell_ids = NULL,
  gene_ids = NULL,
  preprocessing = "none",
  use_second_layer = TRUE,
  nmf_hals_params = params_nmf_hals(),
  nmf_consensus_params = params_nmf_consensus(),
  n_runs = 30L,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `MetaCells` class.

- k:

  Integer. Number of latent factors. Must be at least 2.

- cell_ids:

  Optional character. Cell ids (or meta cell ids) to restrict the NMF
  to. If `NULL`, uses
  [`get_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/get_cells_to_keep.md)
  for `SingleCells` and all meta cells for `MetaCells`.

- gene_ids:

  Optional character. Gene ids to restrict the NMF to. If `NULL`, uses
  [`get_hvg()`](https://gregorlueg.github.io/bixverse/reference/get_hvg.md)
  on the object.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- use_second_layer:

  Boolean. If `TRUE`, runs NMF on the normalised counts (recommended);
  if `FALSE`, on the raw counts.

- nmf_hals_params:

  List, see
  [`params_nmf_hals()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_hals.md).
  The `nmf_init` field is ignored, restarts always use random
  initialisation.

- nmf_consensus_params:

  List, see
  [`params_nmf_consensus()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_consensus.md).

- n_runs:

  Integer. Number of random restarts. Must be at least 2.

- seed:

  Integer. Base random seed. Restart `i` uses `seed + i`, and the
  k-means step is seeded from it too.

- .verbose:

  Boolean or integer. Verbosity.

## Value

A `ConsensusNmfResult` object.

## Details

Use
[`nmf_k_sweep_sc()`](https://gregorlueg.github.io/bixverse/reference/nmf_k_sweep_sc.md)
first if you do not already know `k`.

Memory is the thing to watch. The restarts are dense and all live at
once: roughly `n_cells * k * n_runs` plus `n_runs * k * n_genes` floats
on top of the counts. At 200k cells, `k = 20` and `n_runs = 50` that is
a few hundred megabytes before anything else, which is exactly the
regime where running on `MetaCells` instead is the honest answer.

The density filter is the part that bites. If it leaves fewer than `k`
components the fit errors rather than returning something degenerate.
With few restarts the filter is jumpy, so either raise `n_runs` or set
`density_threshold = 2` in
[`params_nmf_consensus()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_consensus.md)
to switch it off.

## References

Kotliar et al., eLife, 2019
