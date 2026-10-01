# Wrapper function for Bonsai parameters

Parameters for
[`bonsai_sc()`](https://gregorlueg.github.io/bixverse/reference/bonsai_sc.md).
Sanity turns the raw counts into posterior log fold changes with error
bars, Bonsai builds a tree over the cells from those, and the tree is
laid out in 2D for plotting. The Sanity defaults come from sanity-sc-rs,
the search defaults from bonsai-rs.

## Usage

``` r
params_sc_bonsai(
  variance_rule = c("marginalise", "posterior_mean", "max_posterior", "fixed"),
  fixed_variance = NULL,
  variance_bins = 161L,
  variance_min = 0.001,
  variance_max = 50,
  min_signal_to_noise = 1,
  start = c("linkage", "greedy_merge"),
  spr_search = c("approximate", "exact"),
  nni_search = c("approximate", "exact"),
  seed = 42L,
  reroot = TRUE,
  layout = c("equal_angle", "equal_daylight", "dendrogram"),
  hyperbolic = FALSE
)
```

## Arguments

- variance_rule:

  String. How Sanity treats each gene's variance in log fold change.
  `"marginalise"` integrates over the variance grid, `"posterior_mean"`
  and `"max_posterior"` collapse it to one value first (cheaper, less
  accurate error bars), and `"fixed"` uses `fixed_variance` for every
  gene. One of
  `c("marginalise", "posterior_mean", "max_posterior", "fixed")`.
  Defaults to `"marginalise"`.

- fixed_variance:

  Numeric or `NULL`. The variance for `variance_rule = "fixed"`. Ignored
  otherwise. Defaults to `NULL`.

- variance_bins:

  Integer. Bins in the log-spaced variance grid. Defaults to `161L`.

- variance_min:

  Numeric. Smallest variance on the grid. Defaults to `0.001`.

- variance_max:

  Numeric. Largest variance on the grid. Defaults to `50.0`.

- min_signal_to_noise:

  Numeric. Genes whose signal-to-noise ratio falls below this are left
  out of the tree. Defaults to `1.0`.

- start:

  String. Initial topology. `"linkage"` is a Ward linkage over a
  neighbour graph: faster and a better tree on real data.
  `"greedy_merge"` is the paper's star merge, kept for comparisons with
  the published method. One of `c("linkage", "greedy_merge")`. Defaults
  to `"linkage"`.

- spr_search:

  String. Subtree pruning and regrafting. `"approximate"` stays within a
  few nats of `"exact"` on bonsai-rs's benchmarks at up to twice the
  speed. One of `c("approximate", "exact")`. Defaults to
  `"approximate"`.

- nni_search:

  String. Nearest-neighbour interchange. `"approximate"` found the same
  tree as `"exact"` on every bonsai-rs benchmark, and faster. One of
  `c("approximate", "exact")`. Defaults to `"approximate"`.

- seed:

  Integer. Seed for the linkage, SPR and NNI searches. Defaults to
  `42L`.

- reroot:

  Boolean. Reroot the finished tree for display. Changes the drawing,
  not the likelihood. Defaults to `TRUE`.

- layout:

  String. 2D layout. `"equal_daylight"` refines equal angle but falls
  back to it above 2,048 nodes, so it only matters for small trees. Can
  be changed later with
  [`relayout_bonsai()`](https://gregorlueg.github.io/bixverse/reference/relayout_bonsai.md).
  One of `c("equal_angle", "equal_daylight", "dendrogram")`. Defaults to
  `"equal_angle"`.

- hyperbolic:

  Boolean. Project the layout onto the hyperbolic disk. Defaults to
  `FALSE`.

## Value

A named list with the following elements:

- variance_rule - String. How Sanity treats each gene's variance in log
  fold change. `"marginalise"` integrates over the variance grid,
  `"posterior_mean"` and `"max_posterior"` collapse it to one value
  first (cheaper, less accurate error bars), and `"fixed"` uses
  `fixed_variance` for every gene. One of
  `c("marginalise", "posterior_mean", "max_posterior", "fixed")`.
  Defaults to `"marginalise"`.

- fixed_variance - Numeric or `NULL`. The variance for
  `variance_rule = "fixed"`. Ignored otherwise. Defaults to `NULL`.

- variance_bins - Integer. Bins in the log-spaced variance grid.
  Defaults to `161L`.

- variance_min - Numeric. Smallest variance on the grid. Defaults to
  `0.001`.

- variance_max - Numeric. Largest variance on the grid. Defaults to
  `50.0`.

- min_signal_to_noise - Numeric. Genes whose signal-to-noise ratio falls
  below this are left out of the tree. Defaults to `1.0`.

- start - String. Initial topology. `"linkage"` is a Ward linkage over a
  neighbour graph: faster and a better tree on real data.
  `"greedy_merge"` is the paper's star merge, kept for comparisons with
  the published method. One of `c("linkage", "greedy_merge")`. Defaults
  to `"linkage"`.

- spr_search - String. Subtree pruning and regrafting. `"approximate"`
  stays within a few nats of `"exact"` on bonsai-rs's benchmarks at up
  to twice the speed. One of `c("approximate", "exact")`. Defaults to
  `"approximate"`.

- nni_search - String. Nearest-neighbour interchange. `"approximate"`
  found the same tree as `"exact"` on every bonsai-rs benchmark, and
  faster. One of `c("approximate", "exact")`. Defaults to
  `"approximate"`.

- seed - Integer. Seed for the linkage, SPR and NNI searches. Defaults
  to `42L`.

- reroot - Boolean. Reroot the finished tree for display. Changes the
  drawing, not the likelihood. Defaults to `TRUE`.

- layout - String. 2D layout. `"equal_daylight"` refines equal angle but
  falls back to it above 2,048 nodes, so it only matters for small
  trees. Can be changed later with
  [`relayout_bonsai()`](https://gregorlueg.github.io/bixverse/reference/relayout_bonsai.md).
  One of `c("equal_angle", "equal_daylight", "dendrogram")`. Defaults to
  `"equal_angle"`.

- hyperbolic - Boolean. Project the layout onto the hyperbolic disk.
  Defaults to `FALSE`.

## References

de Groot, et al., Nat. Biotechnol., 2026; Breda, et al., Nat.
Biotechnol., 2021.
