# Generate the ligand to target influence matrices

**\[experimental\]** Helper function to generate the ligand to target
influence matrix for the NicheNet like approach.

## Usage

``` r
rs_generate_ligand_target_influence(
  ligand_seeds,
  ppi_network,
  grn_network,
  n_nodes,
  params
)
```

## Arguments

- ligand_seeds:

  List of integer vectors. The 0-indexed seed node(s) per ligand or
  ligand combination.

- ppi_network:

  Named list. Contains the PPI network with the ligand to receptor to
  signalling to TFs. Must contain `from` and `to` (0-indexed node
  indices) and `weight`.

- grn_network:

  Named list. Contains the gene regulatory network with the TF to target
  gene network. Must contain `from` and `to` (0-indexed node indices)
  and `weight`.

- n_nodes:

  Integer. Number of total nodes.

- params:

  Named list. The ligand-target diffusion parameters.

## Value

A dense matrix of ligands x `n_nodes`, rows in `ligand_seeds` order,
with the ligand to target influence scores.
