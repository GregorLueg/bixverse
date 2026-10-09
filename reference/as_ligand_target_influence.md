# Wrap an existing ligand-target matrix as a `LigandTargetInfluence`

Turns a regulatory potential matrix built elsewhere, for example with
nichenetr's `construct_ligand_target_matrix()`, into a
`LigandTargetInfluence` object, so it can be scored with
[`ligand_activity_scores()`](https://gregorlueg.github.io/bixverse/reference/ligand_activity_scores.md).

## Usage

``` r
as_ligand_target_influence(x, ligands_as_cols = FALSE, params = NULL)
```

## Arguments

- x:

  Numeric matrix with unique row and column names.

- ligands_as_cols:

  Boolean. Set to `TRUE` for a genes x ligands matrix, which is the
  nichenetr layout. Defaults to `FALSE` (ligands x genes).

- params:

  Optional list. Parameters used to build `x`, kept for the record.
  Defaults to `NULL`.

## Value

A `LigandTargetInfluence` object. Each ligand is its own seed.

## Examples

``` r
# genes x ligands, as nichenetr returns it
ltm <- matrix(
  c(0.5, 0.0, 0.1, 0.4, 0.3, 0.1),
  nrow = 3,
  dimnames = list(c("G1", "G2", "G3"), c("L1", "L2"))
)
inf <- as_ligand_target_influence(ltm, ligands_as_cols = TRUE)
ligand_activity_scores(inf, gene_sets = list(set_A = "G1"))
#>    gene_set ligand auroc  aupr aupr_corrected   pearson  spearman
#>      <char> <char> <num> <num>          <num>     <num>     <num>
#> 1:    set_A     L1     1     1      0.6666667 0.9819805 0.8660254
#> 2:    set_A     L2     1     1      0.6666667 0.7559289 0.8660254
```
