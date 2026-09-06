# Get the ligand-target influence matrix

Get the ligand-target influence matrix

## Usage

``` r
get_influence(x)

# S3 method for class 'LigandTargetInfluence'
get_influence(x)
```

## Arguments

- x:

  An object holding ligand-target influence results.

## Examples

``` r
# the ligands x genes regulatory potential matrix
ppi <- data.table::data.table(
  from = c("L1", "SIG1", "L2", "SIG2"),
  to = c("SIG1", "TF1", "SIG2", "TF2"),
  weight = 1.0
)
grn <- data.table::data.table(
  from = rep(c("TF1", "TF2"), each = 3),
  to = c("G1", "G2", "G3", "G4", "G5", "G6"),
  weight = 1.0
)
inf <- generate_ligand_target_influence(
  ligand_seeds = list(L1 = "L1", L2 = "L2"),
  ppi_network = ppi,
  grn_network = grn,
  params = params_ligand_target(ltf_cutoff = 0)
)
dim(get_influence(inf))
#> [1]  2 12
```
