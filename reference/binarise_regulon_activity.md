# Binarise regulon activity into on/off calls

The last SCENIC step. Each regulon gets its own threshold derived from
the shape of its AUC distribution across cells, and a cell counts as on
when its score sits strictly above that threshold.

The thresholds come back alongside the calls, so you can inspect them,
override any that look wrong and re-apply the comparison yourself.
SCENIC does the same, it writes them to an editable file between scoring
and assignment.

Regulons flagged as not bimodal fell back to `mean + 2 * sd`. If most of
your regulons land there, the AUC distributions are too flat to
separate, which usually points at the scoring statistic rather than the
thresholding. Check you used `auc_type = "recovery"`, see
[`params_sc_aucell()`](https://gregorlueg.github.io/bixverse/reference/params_sc_aucell.md).

## Usage

``` r
binarise_regulon_activity(
  auc_matrix,
  binarise_params = params_scenic_binarise(),
  .verbose = TRUE
)
```

## Arguments

- auc_matrix:

  Numeric matrix of cells x regulons, or the `ScMatrixRes` returned by
  [`aucell_sc()`](https://gregorlueg.github.io/bixverse/reference/aucell_sc.md).

- binarise_params:

  List. Output of
  [`params_scenic_binarise()`](https://gregorlueg.github.io/bixverse/reference/params_scenic_binarise.md).

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with:

- binary - Logical matrix of cells x regulons, `TRUE` where the regulon
  is on.

- thresholds - data.table with one row per regulon, holding the
  `threshold`, whether it was `bimodal` and the number of cells called
  on.

## References

Aibar, et al., Nat Methods, 2017

## Examples

``` r
# on/off calls for a bimodal and a unimodal regulon
set.seed(7L)
auc <- cbind(
  regulon_a = c(rnorm(50, 0.1, 0.02), rnorm(50, 0.4, 0.02)),
  regulon_b = rnorm(100, 0.2, 0.05)
)
rownames(auc) <- sprintf("cell_%i", 1:100)
binarise_regulon_activity(auc, .verbose = FALSE)$thresholds
#>      regulon threshold bimodal n_cells_on
#>       <char>     <num>  <lgcl>      <num>
#> 1: regulon_a 0.2536428    TRUE         50
#> 2: regulon_b 0.3014215   FALSE          2
```
