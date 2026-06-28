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

A list with predicted_doublets, doublet_scores, threshold,
cluster_labels and detected_doublet_rate.
