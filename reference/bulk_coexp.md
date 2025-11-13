# Bulk RNAseq co-expression modules

Class for applying various co-expression module detection methods on top
of bulk RNAseq data.

## Usage

``` r
bulk_coexp(raw_data, meta_data, variable_info = NULL)
```

## Arguments

- raw_data:

  The raw count matrix. Rows = samples, columns = features.

- meta_data:

  data.table Metadata information on the samples. Expects to have a
  `sample_id` column.

- variable_info:

  data.table. Metadata information on the features. This is an optional
  table.

## Value

Returns the `bulk_coexp` class for further operations.

## Properties

- raw_data:

  A numerical matrix of the provided raw data.

- meta_data:

  A data.table with the meta-information about the samples.

- variable_info:

  An optional data.table containing the variable info.

- processed_data:

  A list in which various types of processed data will be stored.

- outputs:

  A list in which key outputs will be stored.

- params:

  A (nested) list that will store all the parameters of the applied
  function.

- final_results:

  A data.table that will contain the final results.
