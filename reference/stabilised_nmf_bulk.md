# Run stabilised (multi-restart) NMF on a BulkCoExp

Runs `n_runs` HALS-NMF fits with random initialisations and returns the
run with the lowest final reconstruction loss as the primary result,
along with the per-run losses and convergence flags. Useful when the
objective is known to have multiple local minima. Uses `f64` precision.

## Usage

``` r
stabilised_nmf_bulk(
  object,
  k,
  n_runs = 30L,
  preprocessing = c("none", "sd", "sqrt_sd"),
  nmf_hals_params = params_nmf_hals(),
  membership_params = params_module_membership(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

- k:

  Integer. Number of latent factors (modules) per run.

- n_runs:

  Integer. Number of random restarts.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- nmf_hals_params:

  List. Output of
  [`params_nmf_hals()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_hals.md).
  Note that `nmf_init` is ignored for the stabilised variant, which
  always uses random initialisation seeded by `seed + i`.

- membership_params:

  List. Controls how the gene loadings are turned into module
  membership, see
  [`params_module_membership()`](https://gregorlueg.github.io/bixverse/reference/params_module_membership.md).

- seed:

  Integer. Base random seed.

- .verbose:

  Boolean or integer `0L`/`1L`/`2L`. Controls verbosity.

## Value

The class with `final_results` populated from the best run (lowest final
loss) plus stability diagnostics (`losses`, `converged`, `best_idx`,
`w_all_runs`, `h_per_run`).

## Examples

``` r
# ten random restarts, best run kept
syn <- generate_gene_module_data(n_samples = 24L, n_genes = 60L)
# NMF needs a non-negative matrix
mat <- syn$data - min(syn$data)
obj <- BulkCoExp(mat, syn$meta_data)
obj <- preprocess_bulk_coexp(
  obj, hvg = NULL, scaling = FALSE, .verbose = FALSE
)
obj <- stabilised_nmf_bulk(obj, k = 4L, n_runs = 10L, .verbose = FALSE)
get_nmf_stability(obj)$losses
#>  [1] 2.259699 2.258609 2.259214 2.258859 2.258488 2.259210 2.260035 2.261240
#>  [9] 2.260705 2.258274
```
