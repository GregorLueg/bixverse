# Constructor for ligand-target influence results

Stores the ligand x gene influence matrix from the NicheNet-style
computation, the ligand seed groupings, the gene universe, and the
parameters used.

## Usage

``` r
new_ligand_target_influence(
  influence,
  ligand_seeds,
  ligand_names,
  gene_ids,
  params
)
```

## Arguments

- influence:

  Numeric matrix. Rows are ligand seeds, columns are genes.

- ligand_seeds:

  List of character vectors. Each element is a single ligand or a group
  treated as one complex.

- ligand_names:

  Character vector. Row names for the influence matrix.

- gene_ids:

  Character vector. Column names for the influence matrix (full gene
  universe used in the run).

- params:

  List. Parameters used for the run.

## Value

An object of class `LigandTargetInfluence`.
