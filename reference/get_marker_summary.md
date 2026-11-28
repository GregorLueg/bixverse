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
#> 1: cell_type_1 gene_01 0.9880239    0.9531267 0.9520958  0.9531267 0.9541575
#> 2: cell_type_1 gene_02 0.9880239    0.9434446 0.9398312  0.9434446 0.9470580
#> 3: cell_type_1 gene_03 1.0000000    0.9699717 0.9675348  0.9699717 0.9724085
#> 4: cell_type_1 gene_04 0.9580838    0.9179307 0.9086916  0.9179307 0.9271697
#> 5: cell_type_1 gene_05 0.8383234    0.8323234 0.8302915  0.8323234 0.8343554
#> 6: cell_type_1 gene_06 0.9760479    0.9415449 0.9400656  0.9415449 0.9430242
#>    worst_rival min_rank      simes_p    simes_fdr        max_p    max_p_fdr
#>         <char>    <int>        <num>        <num>        <num>        <num>
#> 1: cell_type_3        4 6.581099e-47 8.226374e-46 1.432128e-46 1.790160e-45
#> 2: cell_type_3        5 1.615884e-45 1.615884e-44 3.255928e-44 2.713274e-43
#> 3: cell_type_3        2 1.378589e-50 3.446472e-49 1.178490e-49 2.946226e-48
#> 4: cell_type_2        7 1.081619e-41 7.725847e-41 1.125127e-38 8.036625e-38
#> 5: cell_type_2        9 2.159089e-27 1.199494e-26 5.643505e-27 3.135280e-26
#> 6: cell_type_3        5 9.874756e-45 8.228963e-44 2.109971e-44 2.109971e-43

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
