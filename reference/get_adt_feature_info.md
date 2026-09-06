# Get the ADT feature info

Returns the number of cells expressing a given ADT

## Usage

``` r
get_adt_feature_info(x)

# S3 method for class 'ADTCounts'
get_adt_feature_info(x)
```

## Arguments

- x:

  The object from which to get the ADT feature information

## Value

A data.table with the ADT feature information

## Examples

``` r
# number of cells expressing each protein
adt <- generate_single_cell_test_data_adt()
cell_info <- stats::setNames(
  seq_len(nrow(adt$counts)),
  rownames(adt$counts)
)
adt_clr <- new_adt_counts_clr(adt$counts, cell_info = cell_info)
head(get_adt_feature_info(adt_clr))
#>    feature_idx feature_id   nnz
#>          <int>     <char> <num>
#> 1:           1 protein_01   988
#> 2:           2 protein_02   963
#> 3:           3 protein_03   986
#> 4:           4 protein_04   721
#> 5:           5 protein_05   978
#> 6:           6 protein_06   944
```
