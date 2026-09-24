# Run scDblFinder doublet detection

**\[experimental\]** Implementation of scDblFinder in Rust.

## Usage

``` r
rs_sc_scdblfinder(
  f_path_gene,
  f_path_cell,
  cell_indices,
  params,
  return_features,
  streaming,
  seed,
  verbose
)
```

## Arguments

- f_path_gene:

  String. Path to the gene-based binary file.

- f_path_cell:

  String. Path to the cell-based binary file.

- cell_indices:

  Integer vector (0-indexed).

- params:

  List. scDblFinder parameters from R.

- return_features:

  Boolean. Return the features for the observed cells that are used to
  train the classifier.

- streaming:

  Boolean. Shall the gene data be streamed in for the selection of the
  top genes.

- seed:

  Integer. Seed for reproducibility.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with

- predicted_doublets - Boolean vector (TRUE = doublet).

- doublet_scores - Numerical vector with the classifier probability per
  observed cell.

- cxds_scores - Numerical vector with the cxds scores.

- weighted - Numerical vector with the weighted scores.

- threshold - Threshold used for the doublet calls.

- cluster_labels - Integer vector with the cluster labels from the final
  iteration.

- detected_doublet_rate - Fraction of cells called as doublets.

- selected_genes - Integer vector with the selected gene indices
  (0-indexed!).

- features - If `return_features = TRUE`, a list with `feature_names`,
  `feature_mat` (observed cells x features) and `included_in_training`;
  otherwise an empty list.
