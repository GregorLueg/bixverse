# Bulk RNAseq differential gene expression class

Class for coordinating differential gene expression analyses with
subsequent GSE in a structured format. The filtered counts, library
sizes and normalisation factors are stored in the class;
[`get_dge_list()`](https://gregorlueg.github.io/bixverse/reference/get_dge_list.md)
turns them into an edgeR `DGEList` on demand.

## Usage

``` r
BulkDge(
  raw_counts,
  meta_data,
  variable_info = NULL,
  alternative_gene_id = NULL
)
```

## Arguments

- raw_counts:

  matrix. The raw count matrix. Rows = genes, columns = samples. Note:
  this is different from the
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md)
  class!

- meta_data:

  data.table. Metadata information on the samples. It expects to have a
  column sample_id and case_control column.

- variable_info:

  data.table. Metadata information on the features. This is an optional
  table. Defaults to `NULL`.

- alternative_gene_id:

  String. Optional alternative gene identifier to be used. Must be a
  column of variable_info!

## Value

Returns the `BulkDge` class for further operations.

## Properties

- raw_counts:

  A numerical matrix of the provided raw data.

- meta_data:

  A data.table with the meta-information about the samples.

- variable_info:

  An optional data.table containing the variable info.

- outputs:

  A list in which key outputs will be stored.

- plots:

  A list with the plots that are generated during subsequent QC steps.

- params:

  A (nested) list that will store all the parameters of the applied
  function.

- final_results:

  A list in which final results will be stored.

## Examples

``` r
# DGE class over synthetic bulk counts (genes x samples)
syn <- synthetic_bulk_cor_matrix()
meta <- data.table::data.table(
  sample_id = colnames(syn$counts),
  case_control = rep(c("case", "control"), each = 50)
)
object <- BulkDge(raw_counts = syn$counts, meta_data = meta)
object
#> Bulk differential gene expression class (BulkDge).
#>  Raw counts: 1000 genes x 100 samples.
#>  Meta-data rows: 100.
#>  Variable info provided: FALSE.
#>  Applied steps:
#>   qc_bulk_dge(): FALSE.
#>   normalise_bulk_dge(): FALSE.
#>   batch_correction_bulk_dge(): FALSE.
#>   calculate_pca_bulk_dge(): FALSE.
#>   calculate_dge_limma(): FALSE.
#>   calculate_dge_hedges(): FALSE.
#>   TPM normalisation: FALSE.
#>   FPKM normalisation: FALSE.
```
