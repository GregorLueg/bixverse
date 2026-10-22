# Fit gene trends over Palantir pseudotime

Assigns cells to branches from their fate probabilities, then fits a
landmark Gaussian process with a Matern-5/2 kernel per branch, which is
the Mellon-based estimator of the reference rather than the legacy GAM.
Out comes one smooth curve per gene and branch, on a common pseudotime
grid.

Whatever expression is handed over is what gets fitted. With
`use_magic = FALSE` that is the normalised counts; with `TRUE` it is the
imputed layer
[`run_magic_sc()`](https://gregorlueg.github.io/bixverse/reference/run_magic_sc.md)
wrote. Imputing first smooths twice, which is a defensible presentation
choice and not usually a necessary one.

## Usage

``` r
run_gene_trends_sc(
  object,
  palantir_res,
  features = NULL,
  use_magic = FALSE,
  branch_params = params_sc_branch_selection(),
  gene_trend_params = params_sc_gene_trends(),
  .verbose = TRUE
)
```

## Arguments

- object:

  One of `SingleCells`, `SingleCellsSubset`, `MetaCells` or
  `SingleCellsMultiModal`.

- palantir_res:

  `PalantirRes` class. The output of
  [`run_palantir_sc()`](https://gregorlueg.github.io/bixverse/reference/run_palantir_sc.md),
  run over the same cells.

- features:

  Optional character vector. The genes to fit. If `NULL`, every gene in
  the imputed layer is used, which needs `use_magic = TRUE`.

- use_magic:

  Boolean. Fit the MAGIC imputed layer rather than the normalised
  counts. Needs
  [`run_magic_sc()`](https://gregorlueg.github.io/bixverse/reference/run_magic_sc.md)
  to have been run.

- branch_params:

  List. See
  [`params_sc_branch_selection()`](https://gregorlueg.github.io/bixverse/reference/params_sc_branch_selection.md).

- gene_trend_params:

  List. See
  [`params_sc_gene_trends()`](https://gregorlueg.github.io/bixverse/reference/params_sc_gene_trends.md).

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

A `GeneTrendsRes` S3 object with:

- trends - Long data.table with `branch`, `pseudotime`, `gene` and
  `expression`, ready to facet.

- branch_cells - Named list of the cell names selected for each branch.

- branches - Character vector with the branch names, taken from the
  column names of the fate probability matrix.

- params - List with the parameters the run used.

- run_info - List with `n_cells` and `jitter_used` per branch and the
  grid `resolution`.

## On the defaults

Palantir's pseudotime is min-max scaled to `[0, 1]`, so the reference's
`length_scale` of `1.0` spans the entire domain and its `sigma` of `1.0`
sits at roughly the signal scale of log-normalised expression. The
posterior is prior-dominated: it flattens genuine transient structure
and resolves almost any gene into a smooth monotone or single-peaked
curve. Shorten `length_scale` in
[`params_sc_gene_trends()`](https://gregorlueg.github.io/bixverse/reference/params_sc_gene_trends.md)
before believing a bump.

## References

Setty, et al., Nat. Biotechnol., 2019.

## Examples

``` r
# GP trends for five genes along the Palantir branches
sc <- demo_single_cells()
palantir <- run_palantir_sc(
  sc,
  early_cell = get_knn_obj(sc)$used_cells[1],
  palantir_params = params_sc_palantir(
    knn = 15L,
    num_waypoints = 100L,
    n_eigs = 3L,
    use_early_cell_as_start = TRUE,
    knn_params = list(knn_method = "exhaustive")
  ),
  .verbose = FALSE
)
res <- run_gene_trends_sc(
  sc,
  palantir_res = palantir,
  features = get_gene_names(sc)[1:5],
  gene_trend_params = params_sc_gene_trends(resolution = 50L),
  .verbose = FALSE
)
head(res$trends)
#>      branch pseudotime    gene expression
#>      <fctr>      <num>  <fctr>      <num>
#> 1: cell_099 0.00000000 gene_01   6.259300
#> 2: cell_099 0.02040816 gene_01   6.247292
#> 3: cell_099 0.04081633 gene_01   6.229469
#> 4: cell_099 0.06122449 gene_01   6.205736
#> 5: cell_099 0.08163265 gene_01   6.176020
#> 6: cell_099 0.10204082 gene_01   6.140287

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
