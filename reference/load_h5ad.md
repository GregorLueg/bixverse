# Load in h5ad to `SingleCells`

This function takes an h5ad file and loads the obs and var data into the
DuckDB of the `SingleCells` class and the counts into a Rust-binarised
format for rapid access. During the reading in of the counts, the log
CPM transformation will occur automatically.

## Usage

``` r
load_h5ad(
  object,
  h5_path,
  sc_qc_param = params_sc_min_quality(),
  streaming = 1L,
  raw_count_slot = c("auto", "X", "raw.X", "layers.counts"),
  cell_id_col = NULL,
  batch_size = 1000L,
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- h5_path:

  File path to the h5ad object you wish to load in.

- sc_qc_param:

  List. Output of
  [`params_sc_min_quality()`](https://gregorlueg.github.io/bixverse/reference/params_sc_min_quality.md).
  A list with the following elements:

  - min_unique_genes - Integer. Minimum number of genes to be detected
    in the cell to be included.

  - min_lib_size - Integer. Minimum library size in the cell to be
    included.

  - min_cells - Integer. Minimum number of cells a gene needs to be
    detected to be included.

  - target_size - Float. Target size to normalise to. Defaults to `1e5`.

- streaming:

  Integer. `0L` -\> all cells loaded in memory then transposed (fastest,
  highest memory), `1L` -\> light streaming with cell batching, `2L` -\>
  heavy streaming with memory upper boundaries on the gene side.
  Controls memory pressure during the CSR-to-CSC conversion. Defaults to
  `1L`.

- raw_count_slot:

  Where raw counts live. `"auto"` detects per file via
  [`detect_raw_count_slot()`](https://gregorlueg.github.io/bixverse/reference/detect_raw_count_slot.md);
  otherwise one of `"X"`, `"raw.X"`, `"layers.counts"`.

- cell_id_col:

  Optional string. If a specific column in the h5ad obs data is
  representing the cell identifiers, you can specify it here.

- batch_size:

  Integer. Cell batch size when `streaming = 1L`. Defaults to `1000L`.

- max_genes_in_memory:

  Integer. Maximum genes held in memory at once when `streaming = 2L`.
  Defaults to `2000L`.

- cell_batch_size:

  Integer. Cell batch size when `streaming = 2L`. Defaults to `100000L`.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

It will populate the files on disk and return the class with updated
shape information.
