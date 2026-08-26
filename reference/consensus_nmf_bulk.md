# Run consensus NMF on a BulkCoExp

Runs `n_runs` HALS-NMF restarts, pools their components, drops unstable
ones by local density, k-means clusters the survivors into `k` groups
and refits the partner factor against the per-cluster median. This is
cNMF: the answer is the structure the restarts agree on, and the mean
silhouette of those clusters (`stability`) tells you how much they
agreed.

Prefer this over
[`stabilised_nmf_bulk()`](https://gregorlueg.github.io/bixverse/reference/stabilised_nmf_bulk.md),
which just picks the lowest-loss restart and tells you nothing about
whether that restart is reproducible.

## Usage

``` r
consensus_nmf_bulk(
  object,
  k,
  n_runs = 30L,
  preprocessing = c("none", "sd", "sqrt_sd"),
  nmf_hals_params = params_nmf_hals(),
  nmf_consensus_params = params_nmf_consensus(),
  membership_params = params_module_membership(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).
  Ideally, you should run
  [`preprocess_bulk_coexp()`](https://gregorlueg.github.io/bixverse/reference/preprocess_bulk_coexp.md)
  before applying this function. The NMF requires non-negative input, so
  the pre-processing should either skip scaling or use a strictly
  non-negative representation.

- k:

  Integer. Number of latent factors (modules). Must be at least 2.

- n_runs:

  Integer. Number of random restarts. Must be at least 2.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`. Applied inside the Rust
  kernel to the input matrix before the HALS updates.

- nmf_hals_params:

  List. Output of
  [`params_nmf_hals()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_hals.md).
  The `nmf_init` field is ignored, restarts always use random
  initialisation.

- nmf_consensus_params:

  List. Output of
  [`params_nmf_consensus()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_consensus.md).
  Controls the clustering of the pooled components.

- membership_params:

  List. Controls how the gene loadings are turned into module
  membership, see
  [`params_module_membership()`](https://gregorlueg.github.io/bixverse/reference/params_module_membership.md).

- seed:

  Integer. Base random seed. Restart `i` uses `seed + i`, and the
  k-means step is seeded from it too.

- .verbose:

  Boolean or integer `0L`/`1L`/`2L`. Controls verbosity.

## Value

The class with `final_results` populated from the consensus fit, with
the clustering diagnostics under `diagnostics` and the fit parameters
under `params$nmf_fit`.

## Details

Use
[`nmf_k_sweep_bulk()`](https://gregorlueg.github.io/bixverse/reference/nmf_k_sweep_bulk.md)
first if you do not already know `k`.

The density filter is the part that bites. Components whose mean cosine
distance to their neighbours exceeds `density_threshold` are dropped,
and if that leaves fewer than `k` survivors the fit errors rather than
returning something degenerate. With few restarts the filter is jumpy,
so either raise `n_runs` or set `density_threshold = 2` to switch it
off.

## References

Kotliar et al., eLife, 2019
