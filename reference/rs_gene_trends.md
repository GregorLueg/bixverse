# Fit Palantir gene trends over pseudotime

**\[experimental\]** Selects the cells belonging to each branch from the
fate probabilities, then fits a landmark Gaussian process with a
Matern-5/2 kernel per branch. This is the Mellon-based estimator of the
reference, not the legacy GAM.

Whatever expression matrix is handed over is what gets fitted. Nothing
here imputes; that decision belongs to the caller.

The defaults are prior-dominated. Palantir's pseudotime is min-max
scaled to `[0, 1]`, so a `length_scale` of `1.0` spans the whole domain
and a `sigma` of `1.0` sits at roughly the signal scale of
log-normalised expression. That resolves almost any gene into a smooth
monotone or single-peaked curve. Shorten `length_scale` before believing
a bump.

## Usage

``` r
rs_gene_trends(
  expression,
  pseudotime,
  branch_probs,
  branch_params,
  gene_trend_params
)
```

## Arguments

- expression:

  Numerical matrix of cells x genes. All cells, in the same row order as
  `pseudotime`, not just the branch members.

- pseudotime:

  Numerical vector. Pseudotime per cell.

- branch_probs:

  Numerical matrix of cells x fates with the fate probabilities. Rows
  need not sum to one.

- branch_params:

  List. Parameter list, see
  [`params_sc_branch_selection()`](https://gregorlueg.github.io/bixverse/reference/params_sc_branch_selection.md).

- gene_trend_params:

  List. Parameter list, see
  [`params_sc_gene_trends()`](https://gregorlueg.github.io/bixverse/reference/params_sc_gene_trends.md).

## Value

A list with

- trends - List of numerical matrices, one per branch, each of
  resolution x genes.

- grids - List of numerical vectors with the pseudotime grid per branch,
  running from the branch minimum to its maximum.

- branch_cells - List of integer vectors with the cell indices
  (0-indexed!) selected for each branch.

- n_cells - Integer vector with the cell count per branch.

- jitter_used - Numerical vector with the jitter each branch's Cholesky
  needed.

## References

Setty, et al., Nat. Biotechnol., 2019.
