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

  List. Contains the indices of the seeds, i.e., ligands.

- ppi_network:

  Named list. Contains the PPI network with the ligand to receptor to
  signalling to TFs. Must contain from (indices), to (indices), and edge
  weights.

- grn_network:

  Named list. Contains the gene regulatory network with the TF to target
  gene network. Must contain from (indices), to (indices), and edge
  weights.

- n_nodes:

  Integer. Number of total nodes.

- params:

  Named list.

## Value

A dense matrix of ligands x genes that contains the influence scores of
each
