# Pre-scan multiple 10x CellRanger h5 files for multi-sample loading

Walks each input file, reads the per-file feature ids, restricts to the
chosen `feature_type` (V3 only; ignored for V2), builds the intersection
or union universe of gene ids, and returns the file tasks expected by
[`load_multi_tenx_h5()`](https://gregorlueg.github.io/bixverse/reference/load_multi_tenx_h5.md).

V2 and V3 files can be mixed in a single batch. Gene matching is done on
the feature `id` (ensembl-style) since gene symbols collide.

## Usage

``` r
prescan_tenx_h5_files(
  h5_paths,
  feature_type = "Gene Expression",
  gene_universe = c("intersection", "union"),
  .verbose = TRUE
)
```

## Arguments

- h5_paths:

  Character vector of file paths to 10x h5 files. If names are provided,
  these will be used as experimental identifiers.

- feature_type:

  String. Modality to keep across all files (V3 only; ignored for V2).
  Defaults to `"Gene Expression"`.

- gene_universe:

  One of `"intersection"` or `"union"`.

- .verbose:

  Boolean. Controls verbosity.

## Value

A list with:

- universe - Character vector of gene ids in the universe.

- universe_size - Length of the universe.

- file_tasks - Named list of per-file task structures, each containing
  `exp_id`, `h5_path`, `version`, `no_cells`, `no_genes`, `feature_type`
  and `gene_local_to_universe` (integer vector, `NA` for features
  outside the universe / non-target modality, 0-indexed).
