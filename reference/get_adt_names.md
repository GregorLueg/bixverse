# Get the ADT feature names

Get the main ADT feature names

## Usage

``` r
get_adt_names(x)

# S3 method for class 'ADTCounts'
get_adt_names(x)
```

## Arguments

- x:

  An object to get the gene names from.

## Value

The primary ADT feature identifiers stored in the class.

## Examples

``` r
# protein identifiers held by an `ADTCounts`
adt <- generate_single_cell_test_data_adt()
cell_info <- stats::setNames(
  seq_len(nrow(adt$counts)),
  rownames(adt$counts)
)
adt_clr <- new_adt_counts_clr(adt$counts, cell_info = cell_info)
head(get_adt_names(adt_clr))
#> [1] "protein_01" "protein_02" "protein_03" "protein_04" "protein_05"
#> [6] "protein_06"
```
