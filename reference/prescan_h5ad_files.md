# Pre-scan multiple h5ad files for multi-sample loading

Pre-scan multiple h5ad files for multi-sample loading

## Usage

``` r
prescan_h5ad_files(
  h5_paths,
  gene_universe = c("intersection", "union"),
  var_index = "_index",
  raw_count_slot = c("auto", "X", "raw.X", "layers.counts"),
  .verbose = TRUE
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

- raw_count_slot:

  Where raw counts live. `"auto"` detects per file via
  [`detect_raw_count_slot()`](https://gregorlueg.github.io/bixverse/reference/detect_raw_count_slot.md);
  otherwise one of `"X"`, `"raw.X"`, `"layers.counts"`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with:

- universe - Character vector of gene names in the universe

- file_tasks - List of per-file task structures, each containing:
  exp_id, h5_path, cs_type, no_cells, no_genes, gene_local_to_universe
  (integer vector, NA for genes not in universe, 0-indexed into
  universe)

## Examples

``` r
# build the gene universe across two files before a multi-file load
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
)
files <- c(a = tempfile(fileext = ".h5ad"), b = tempfile(fileext = ".h5ad"))
for (f in files) {
  write_h5ad_sc(f, data$counts, data$obs, data$var, .verbose = FALSE)
}
tasks <- prescan_h5ad_files(h5_paths = files, .verbose = FALSE)
tasks$universe_size
#> [1] 40

unlink(files)
```
