# Wrapper function for parameters for NEBULA

Parameters for the NEBULA negative binomial gamma mixed model,
implemented in Rust via the `edge-rs` crate and ported from the `nebula`
package's own C++. Defaults are the R package's own.

NEBULA splits the variance into a subject-level random effect and a
cell-level overdispersion. Run it on meta cells and the cell-level term
becomes the spread between aggregates within a subject rather than
between cells, so it is smaller and absorbs whatever the aggregation
smoothed away. The subject-level term keeps its meaning either way.

## Usage

``` r
params_nebula(
  nebula_method = c("ln", "hl"),
  min_sigma = 1e-04,
  min_phi = 1e-04,
  max_sigma = 10,
  max_phi = 1000,
  cutoff_cell = 20,
  kappa = 800,
  cpc = 0.005,
  mincp = 5L,
  reml = FALSE,
  eps = 1e-06,
  gene_batch_size = 1000L,
  shrink_dispersion = TRUE
)
```

## Arguments

- nebula_method:

  String. Which variant to run. One of `c("ln", "hl")`. Defaults to
  `"ln"`. NEBULA downgrades `"ln"` to `"hl"` below 30 cells per subject,
  as the R package does.

- min_sigma:

  Numeric. Lower bound on the subject-level overdispersion. Defaults to
  `1e-4`.

- min_phi:

  Numeric. Lower bound on the cell-level overdispersion. Defaults to
  `1e-4`.

- max_sigma:

  Numeric. Upper bound on the subject-level overdispersion. Defaults to
  `10`.

- max_phi:

  Numeric. Upper bound on the cell-level overdispersion. Defaults to
  `1000`.

- cutoff_cell:

  Numeric. Refit both overdispersions when the product of the cells per
  subject and the estimated `phi` falls below this. Defaults to `20`.

- kappa:

  Numeric. Threshold on NEBULA's `kappa_obs` above which the
  subject-level overdispersion from stage one is trusted as is. Defaults
  to `800`.

- cpc:

  Numeric. Drop a gene whose mean count per cell is at most this.
  Defaults to `0.005`.

- mincp:

  Integer. Drop a gene expressed in fewer than this many cells. Defaults
  to `5L`.

- reml:

  Boolean. Estimate the overdispersions by restricted maximum
  likelihood. Defaults to `FALSE`. The R package only honours this for
  `NBLMM`, which the Rust port does not implement, so this arm has not
  been validated against an R reference. Leave it off unless you know
  why you want it.

- eps:

  Numeric. Absolute stopping tolerance for the optimiser. Defaults to
  `1e-6`.

- gene_batch_size:

  Integer. Genes read and fitted per batch. Bounds how much of the store
  is resident at once and changes nothing about the answer, since NEBULA
  is gene-independent. Defaults to `1000L`.

- shrink_dispersion:

  Boolean. Shrink the cell-level overdispersions towards an empirical
  Bayes prior once the sweep is done. Defaults to `TRUE`.

## Value

A list with the NEBULA parameters.

## References

He, et al., Commun Biol, 2021
