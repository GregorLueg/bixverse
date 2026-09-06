# Get the per-gene marker summaries across all rivals

Returns the summaries a marker is judged on: the AUROC of the reference
against its rivals reduced to one row per gene and reference group. Rank
on `median_auroc` for a marker that survives a single closely related
rival, or on `min_auroc` when it has to beat every rival unambiguously.

## Usage

``` r
get_marker_summary(x)

# S3 method for class 'ScSpecificMarkers'
get_marker_summary(x)
```

## Arguments

- x:

  `ScSpecificMarkers` object.

## Value

A copy of the summary data.table.

## Examples

``` r
# per gene summaries across every rival group
sc <- demo_single_cells()
res <- find_specific_markers_sc(
  sc,
  column_of_interest = "cell_grp",
  .verbose = FALSE
)
head(get_marker_summary(res))
#>        ref_grp gene_id  prop_ref median_auroc min_auroc mean_auroc max_auroc
#>         <char>  <char>     <num>        <num>     <num>      <num>     <num>
#> 1: cell_type_1 gene_01 0.9880239    0.9531087 0.9520777  0.9531087 0.9541396
#> 2: cell_type_1 gene_02 0.9880239    0.9434354 0.9397771  0.9434354 0.9470938
#> 3: cell_type_1 gene_03 1.0000000    0.9699715 0.9674987  0.9699715 0.9724444
#> 4: cell_type_1 gene_04 0.9580838    0.9179306 0.9087095  0.9179306 0.9271517
#> 5: cell_type_1 gene_05 0.8383234    0.8323145 0.8302556  0.8323145 0.8343734
#> 6: cell_type_1 gene_06 0.9760479    0.9415538 0.9400476  0.9415538 0.9430600
#>    worst_rival min_rank      simes_p    simes_fdr        max_p    max_p_fdr
#>         <char>    <int>        <num>        <num>        <num>        <num>
#> 1: cell_type_3        4 6.632426e-47 8.290532e-46 1.443307e-46 1.804134e-45
#> 2: cell_type_3        5 1.589201e-45 1.589201e-44 3.333359e-44 2.777799e-43
#> 3: cell_type_3        2 1.354582e-50 3.386456e-49 1.197853e-49 2.994632e-48
#> 4: cell_type_2        7 1.089790e-41 7.784218e-41 1.116568e-38 7.975489e-38
#> 5: cell_type_2        9 2.145260e-27 1.191811e-26 5.714230e-27 3.174572e-26
#> 6: cell_type_3        5 9.715047e-45 8.095872e-44 2.126360e-44 2.126360e-43

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
