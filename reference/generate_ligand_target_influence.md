# Generate the ligand to target influence matrix

Computes the NicheNet-style ligand to target gene influence matrix.
Builds the gene universe as the union of all symbols across both
networks, remaps to 0-indexed integer node IDs for the Rust side, and
wraps the resulting matrix in a `LigandTargetInfluence` object.

## Usage

``` r
generate_ligand_target_influence(
  ligand_seeds,
  ppi_network,
  grn_network,
  params = params_ligand_target()
)
```

## Arguments

- ligand_seeds:

  List of character vectors. Each element is one ligand symbol or a
  group treated as one complex (e.g. `list(c("TGFB1", "TGFB2"))`).
  Optionally named; if unnamed, row names default to the symbols joined
  with `+`.

- ppi_network:

  data.table with columns `from`, `to`, `weight` (character, character,
  numeric). Protein-protein / signalling layer.

- grn_network:

  data.table with columns `from`, `to`, `weight` (character, character,
  numeric). Gene regulatory layer.

- params:

  List. As returned by
  [`params_ligand_target()`](https://gregorlueg.github.io/bixverse/reference/params_ligand_target.md).

## Value

A `LigandTargetInfluence` object.

## Examples

``` r
# two disjoint signalling components over a toy network
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
inf
#> LigandTargetInfluence
#>   No ligand seeds:    2
#>   No genes:           12
#>   Damping factor:     0.500
#>   Max iter:           1000
#>   Secondary targets:  FALSE
```
