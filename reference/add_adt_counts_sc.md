# Add ADT counts to `SingleCellsMultiModal`

This method allows you to add ADT counts to a `SingleCellsMultiModal`.
Assumes the transcriptomics data has already been ingested and the cells
to keep are known. The method subsets the ADT counts to the kept cells,
applies the requested normalisation (CLR or DSB), populates the
`"var_adt"` table in the DuckDB, and attaches an `ADTCounts` to the
class.

## Usage

``` r
add_adt_counts_sc(object, adt_counts, method = c("clr", "dsb"), ...)
```

## Arguments

- object:

  `SingleCellsMultiModal` class.

- adt_counts:

  Numeric matrix. Cells x features matrix of raw ADT counts.

- method:

  String. One of `c("clr", "dsb")`. Normalisation method.

- ...:

  Additional arguments forwarded to the normalisation constructor. For
  `method = "clr"`: `seurat_clr`, `clean_clr_counts`, `percentile`. For
  `method = "dsb"`: `empty_drops`, `isotype_names`, `dsb_params`,
  `scale_factor`, `seed`, `verbose`. See
  [`new_adt_counts_clr()`](https://gregorlueg.github.io/bixverse/reference/new_adt_counts_clr.md)
  and
  [`new_adt_counts_dsb()`](https://gregorlueg.github.io/bixverse/reference/new_adt_counts_dsb.md).

## Value

Returns a `SingleCellsMultiModal` with the ADT data added.

## Examples

``` r
# a CLR-normalised ADT layer on top of an ingested RNA modality
rna <- generate_single_cell_test_data()
adt <- generate_single_cell_test_data_adt()
dir <- tempfile("bixverse_mm")
dir.create(dir)
object <- load_r_data(
  SingleCellsMultiModal(dir_data = dir),
  counts = rna$counts,
  obs = rna$obs,
  var = rna$var,
  sc_qc_param = params_sc_min_quality(min_unique_genes = 5L),
  .verbose = FALSE
)
object <- add_adt_counts_sc(object, adt_counts = adt$counts, method = "clr")
get_adt_names(object)
#>  [1] "protein_01" "protein_02" "protein_03" "protein_04" "protein_05"
#>  [6] "protein_06" "protein_07" "protein_08" "protein_09" "protein_10"
#> [11] "protein_11" "protein_12" "protein_13" "protein_14" "protein_15"

unlink(dir, recursive = TRUE, force = TRUE)
```
