# Load in data from a `SingleCellExperiment`

Brings a Bioconductor `SingleCellExperiment` into a `SingleCells`
object. `colData` becomes the obs table, `rowData` becomes the var
table, and the chosen assay goes through the same Rust quality control
and normalisation every other loader uses.

The assay has to hold raw counts. Plenty of objects in the wild ship
only `logcounts`, and a negative binomial cannot model those, so pick
the right one rather than letting the default find whatever is there.

`reducedDims` and `altExps` are not carried over. Run the embedding on
this side, and use
[`SingleCellsMultiModal()`](https://gregorlueg.github.io/bixverse/reference/SingleCellsMultiModal.md)
for ADT.

## Usage

``` r
load_sce(
  object,
  sce,
  assay_name = "counts",
  sc_qc_param = params_sc_min_quality(),
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

- sce:

  `SingleCellExperiment` class you want to transform.

- assay_name:

  String. Which assay holds the raw counts. Defaults to `"counts"`.

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

  Integer. CSR-to-CSC conversion mode. `0L` -\> in-memory (fastest,
  highest memory), `1L` -\> light streaming with cell batching, `2L` -\>
  heavy streaming with memory upper boundaries. Defaults to `1L`.

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
