# Run non-negative matrix factorisation on a BulkCoExp

Fits a single HALS-NMF `V ~ W H` to the bulk expression matrix with a
fixed number of components `k`. `V` is samples x features (matches the
layout of `raw_data` / `processed_data` on
[`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md)).
The result is rearranged so that:

- `gene_loadings` (features x k) captures per-module gene contributions.

- `sample_activity` (samples x k) captures per-module sample activity.

Module membership is derived by keeping the upper tail of each
component's gene loadings, see
[`params_module_membership()`](https://gregorlueg.github.io/bixverse/reference/params_module_membership.md).
A gene can belong to several modules, and a gene in no tail belongs to
none.

## Usage

``` r
nmf_bulk(
  object,
  k,
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
  Ideally, you should run
  [`preprocess_bulk_coexp()`](https://gregorlueg.github.io/bixverse/reference/preprocess_bulk_coexp.md)
  before applying this function. The NMF requires non-negative input, so
  the pre-processing should either skip scaling or use a strictly
  non-negative representation.

- k:

  Integer. Number of latent factors (modules).

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`. Applied inside the Rust
  kernel to the input matrix before the HALS updates.

- nmf_hals_params:

  List. Output of
  [`params_nmf_hals()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_hals.md):

  - max_iter - Integer. Maximum number of HALS iterations.

  - tol - Float. Convergence tolerance.

  - eps - Float. Numerical floor.

  - check_every - Integer. Convergence check interval.

  - nmf_init - String. One of `c("nndsvd", "svd", "random")`.

- membership_params:

  List. Controls how the gene loadings are turned into module
  membership, see
  [`params_module_membership()`](https://gregorlueg.github.io/bixverse/reference/params_module_membership.md).
  Membership is not exclusive: a gene loading strongly on several
  components appears in several modules, and a gene in no tail appears
  in none.

- seed:

  Integer. Random seed for the NMF initialisation. Defaults to `42L`.

- .verbose:

  Boolean or integer `0L`/`1L`/`2L`. Controls verbosity.

## Value

The class with `final_results` populated (see description) and the fit
parameters stored under `params$nmf_fit`, plus
`params$detection_method = "nmf-based"`.

## References

Cichocki & Phan, IEICE Trans., 2009.
