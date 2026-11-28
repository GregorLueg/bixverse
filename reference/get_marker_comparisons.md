# Get the per-rival marker statistics

Returns the statistics that the summaries of
[`get_marker_summary()`](https://gregorlueg.github.io/bixverse/reference/get_marker_summary.md)
are built from, i.e., one row per gene, reference group and rival.
Useful to find out which rival a gene fails against.

## Usage

``` r
get_marker_comparisons(x)

# S3 method for class 'ScSpecificMarkers'
get_marker_comparisons(x)
```

## Arguments

- x:

  `ScSpecificMarkers` object.

## Value

A copy of the per comparison data.table.

## Examples

``` r
# the per rival statistics the summaries are built from
sc <- demo_single_cells()
res <- find_specific_markers_sc(
  sc,
  column_of_interest = "cell_grp",
  .verbose = FALSE
)
head(get_marker_comparisons(res))
#>        ref_grp   rival_grp gene_id     auroc      lfc  prop_ref prop_rival
#>         <char>      <char>  <char>     <num>    <num>     <num>      <num>
#> 1: cell_type_1 cell_type_2 gene_01 0.9541575 3.339092 0.9880239  0.6946108
#> 2: cell_type_1 cell_type_2 gene_02 0.9470580 3.292497 0.9880239  0.6886228
#> 3: cell_type_1 cell_type_2 gene_03 0.9724085 3.461720 1.0000000  0.6946108
#> 4: cell_type_1 cell_type_2 gene_04 0.9086916 2.915461 0.9580838  0.6886228
#> 5: cell_type_1 cell_type_2 gene_05 0.8302915 2.709359 0.8383234  0.4431138
#> 6: cell_type_1 cell_type_2 gene_06 0.9430242 3.113653 0.9760479  0.7005988
#>    z_scores     p_values          fdr
#>       <num>        <num>        <num>
#> 1: 14.38336 3.290550e-47 4.113187e-46
#> 2: 14.16017 8.079419e-46 8.079419e-45
#> 3: 14.95812 6.892943e-51 1.723236e-49
#> 4: 12.95331 1.125127e-38 8.036625e-38
#> 5: 10.69039 5.643505e-27 3.135280e-26
#> 6: 14.03240 4.937378e-45 4.114482e-44

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
