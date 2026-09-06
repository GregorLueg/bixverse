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
#> 1: cell_type_1 cell_type_2 gene_01 0.9541396 3.339259 0.9880239  0.6946108
#> 2: cell_type_1 cell_type_2 gene_02 0.9470938 3.292560 0.9880239  0.6886228
#> 3: cell_type_1 cell_type_2 gene_03 0.9724444 3.461920 1.0000000  0.6946108
#> 4: cell_type_1 cell_type_2 gene_04 0.9087095 2.915536 0.9580838  0.6886228
#> 5: cell_type_1 cell_type_2 gene_05 0.8302556 2.709311 0.8383234  0.4431138
#> 6: cell_type_1 cell_type_2 gene_06 0.9430600 3.113609 0.9760479  0.7005988
#>    z_scores     p_values          fdr
#>       <num>        <num>        <num>
#> 1: 14.38283 3.316213e-47 4.145266e-46
#> 2: 14.16134 7.946003e-46 7.946003e-45
#> 3: 14.95929 6.772912e-51 1.693228e-49
#> 4: 12.95389 1.116568e-38 7.975489e-38
#> 5: 10.68923 5.714230e-27 3.174572e-26
#> 6: 14.03355 4.857523e-45 4.047936e-44

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
