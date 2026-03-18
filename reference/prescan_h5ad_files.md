# Pre-scan multiple h5ad files for multi-sample loading

Pre-scan multiple h5ad files for multi-sample loading

## Usage

``` r
prescan_h5ad_files(
  h5_paths,
  gene_universe = c("intersection", "union"),
  var_index = "_index"
)
```

## Arguments

- h5_paths:

  Character vector of file paths to h5ad files. If names are provided,
  these will be used as experimental identifiers.

- gene_universe:

  One of "intersection" or "union".

- var_index:

  String. The name within the h5ad var part in which the variable names
  are stored. Defaults to `"_index"`.

## Value

A list with:

- universe - Character vector of gene names in the universe

- file_tasks - List of per-file task structures, each containing:
  exp_id, h5_path, cs_type, no_cells, no_genes, gene_local_to_universe
  (integer vector, NA for genes not in universe, 0-indexed into
  universe)
