# Change the primary gene identifier of bulk_dge

Changes the primary gene identifier in the bulk_dge class. To do so, you
need to either provide a `variable_info` data.table with the alternative
gene identifier you wish to use or it exists already in the object
itself. If it exists in the object, that variable_info will be used.

## Usage

``` r
change_gene_identifier(object, alternative_gene_id, variable_info = NULL)
```

## Arguments

- object:

  `bulk_dge` class, see [bulk_dge](bulk_dge.md).

- alternative_gene_id:

  String. The column containing the alternative gene identifier. Must be
  present in the provided `variable_info` data.table or within the class
  attributes.

- variable_info:

  Optional data.table with variable information. If `variable_info` is
  in an attribute of the class, that one will be used.

## Value

The class with modified primary gene identifier.
