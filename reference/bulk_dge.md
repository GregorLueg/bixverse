# Bulk RNAseq differential gene expression class

Class for coordinating differential gene expression analyses with
subsequent GSE in a structured format. Additionally, the class will
store the counts in
[`edgeR::DGEList()`](https://rdrr.io/pkg/edgeR/man/DGEList.html) for
subsequent processing.

## Usage

``` r
bulk_dge(
  raw_counts,
  meta_data,
  variable_info = NULL,
  alternative_gene_id = NULL
)
```

## Arguments

- raw_counts:

  matrix. The raw count matrix. Rows = genes, columns = samples. Note:
  this is different from the [`bulk_coexp()`](bulk_coexp.md) class!

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

Returns the `bulk_coexp` class for further operations.

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
