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
