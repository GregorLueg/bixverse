# Stream in h5ad to `single_cell_exp`

This function takes an h5ad file and loads (via streaming) the obs and
var data into the DuckDB of the `single_cell_exp` class and the counts
into a Rust-binarised format for rapid access. During the reading in of
the counts, the log CPM transformation will occur automatically. This
function is specifically designed to deal with larger amounts of data
and is slower than [`load_h5ad()`](load_h5ad.md).

## Usage

``` r
stream_h5ad(
  object,
  h5_path,
  sc_qc_param = params_sc_min_quality(),
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- h5_path:

  File path to the h5ad object you wish to load in.

- sc_qc_param:

  List. Output of [`params_sc_min_quality()`](params_sc_min_quality.md).
  A list with the following elements:

  - min_unique_genes - Integer. Minimum number of genes to be detected
    in the cell to be included.

  - min_lib_size - Integer. Minimum library size in the cell to be
    included.

  - min_cells - Integer. Minimum number of cells a gene needs to be
    detected to be included.

  - target_size - Float. Target size to normalise to. Defaults to `1e5`.

- max_genes_in_memory:

  Integer. How many genes shall be held in memory at a given point.
  Defaults to `2000L`.

- cell_batch_size:

  Integer. How big are the batch sizes for the cells in the
  transformation from the cell-based to gene-based format. Defaults to
  `100000L`.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

It will populate the files on disk and return the class with updated
shape information.
