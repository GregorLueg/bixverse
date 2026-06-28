# Load in a 10x CellRanger h5 file to `SingleCells`

Loads the gene-expression modality from a CellRanger v2/v3 h5 file. The
counts go into the Rust-binarised format (with log normalisation applied
on read) and the barcodes/features into the DuckDB. Non-gene modalities
(e.g. Antibody Capture) are filtered out via `feature_type`.

## Usage

``` r
load_tenx_h5(
  object,
  h5_path,
  sc_qc_param = params_sc_min_quality(),
  feature_type = "Gene Expression",
  streaming = 1L,
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

  File path to the 10x h5 file.

- sc_qc_param:

  List. Output of
  [`params_sc_min_quality()`](https://gregorlueg.github.io/bixverse/reference/params_sc_min_quality.md).

- feature_type:

  String. Modality to keep. Defaults to `"Gene Expression"`. Ignored for
  v2 (single modality).

- streaming:

  Integer. CSR-to-CSC conversion mode. `0L` -\> in-memory, `1L` -\>
  light streaming, `2L` -\> heavy streaming. Defaults to `1L`.

- batch_size:

  Integer. Cell batch size when `streaming = 1L`. Defaults to `1000L`.

- max_genes_in_memory:

  Integer. Maximum genes held in memory at once when `streaming = 2L`.
  Defaults to `2000L`.

- cell_batch_size:

  Integer. Cell batch size when `streaming = 2L`. Defaults to `100000L`.

- .verbose:

  Boolean.

## Value

The class with updated shape information.
