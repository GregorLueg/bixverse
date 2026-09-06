# Change the primary gene identifier of BulkDge

Changes the primary gene identifier in the BulkDge class. To do so, you
need to either provide a `variable_info` data.table with the alternative
gene identifier you wish to use or it exists already in the object
itself. If it exists in the object, that variable_info will be used.

## Usage

``` r
change_gene_identifier(object, alternative_gene_id, variable_info = NULL)
```

## Arguments

- object:

  `BulkDge` class, see
  [`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md).

- alternative_gene_id:

  String. The column containing the alternative gene identifier. Must be
  present in the provided `variable_info` data.table or within the class
  attributes.

- variable_info:

  Optional data.table with variable information. If `variable_info` is
  in an attribute of the class, that one will be used.

## Value

The class with modified primary gene identifier.

## Examples

``` r
# swap the row names over to a symbol column of the variable info
syn <- synthetic_bulk_cor_matrix()
meta <- data.table::data.table(sample_id = colnames(syn$counts))
var_info <- data.table::data.table(
  var_id = rownames(syn$counts),
  symbol = sprintf("SYM%i", seq_len(nrow(syn$counts)))
)
object <- BulkDge(raw_counts = syn$counts, meta_data = meta)
object <- change_gene_identifier(object, "symbol", var_info)
head(rownames(S7::prop(object, "raw_counts")))
#> [1] "SYM1" "SYM2" "SYM3" "SYM4" "SYM5" "SYM6"
```
