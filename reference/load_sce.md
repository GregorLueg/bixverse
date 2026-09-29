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
  csc_mem_gb = NULL,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
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

- csc_mem_gb:

  Optional numeric. Memory in GB for the buffers of the cell-to-gene
  (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL` (default)
  converts in one pass and holds the whole matrix. Set a cap for large
  data sets; every extra phase re-reads the cell file once.

- streaming, batch_size, max_genes_in_memory, cell_batch_size:

  Replaced by `csc_mem_gb` and ignored. **\[deprecated\]**

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

It will populate the files on disk and return the class with updated
shape information.

## Examples

``` r
# \donttest{
# colData becomes obs, rowData becomes var, the counts assay gets normalised
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
)
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = as(Matrix::t(data$counts), "CsparseMatrix")),
  colData = data.frame(data$obs, row.names = data$obs$cell_id),
  rowData = data.frame(data$var, row.names = data$var$gene_id)
)
dir_data <- tempfile("sc_sce")
dir.create(dir_data, recursive = TRUE)
sc <- load_sce(
  object = SingleCells(dir_data = dir_data),
  sce = sce,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 5L,
    min_lib_size = 25L,
    min_cells = 5L
  ),
  .verbose = FALSE
)
dim(sc)
#> [1] 200  40

unlink(dir_data, recursive = TRUE, force = TRUE)
# }
```
