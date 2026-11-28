# Count fragments over gene body plus promoter into a gene activity store

Signac's `GeneActivity()`: per gene, the longest span of its biotype,
under `max_width`, extended `upstream` of the TSS; each overlapping
fragment counts once. The result is a
[`SingleCells()`](https://gregorlueg.github.io/bixverse/reference/SingleCells.md)
over every cell of `object`, in its order and with its obs and cell
filter, so it serves as the query of
[`transfer_labels_seurat_sc()`](https://gregorlueg.github.io/bixverse/reference/transfer_labels_seurat_sc.md)
and of the RNA tools. Methods that assume RNA counts refuse it or warn.
Written to the default directory, later cell filters on `object` reach
it; reopen it with
[`get_gene_activity()`](https://gregorlueg.github.io/bixverse/reference/get_gene_activity.md).

## Usage

``` r
gene_activity_sc(
  object,
  fragments,
  gtf,
  dir_out = NULL,
  ga_params = params_atac_gene_activity(),
  target_size = NULL,
  barcode_map = NULL,
  overwrite = FALSE,
  max_mem_gb = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class whose cell ids are the fragment barcodes (a
  `SingleCellsMultiModal` or a `SingleCellsAtac`).

- fragments:

  String. The fragment file (BGZF, gzip or plain).

- gtf:

  String. Gene annotation (GENCODE, Ensembl, or a transcript-only GTF).

- dir_out:

  Optional string. Output directory. `NULL` uses `gene_activity` inside
  the experiment directory.

- ga_params:

  List. Output of
  [`params_atac_gene_activity()`](https://gregorlueg.github.io/bixverse/reference/params_atac_gene_activity.md).

- target_size:

  Optional numeric. Target of the log-CPM norm. `NULL` uses the median
  RNA `lib_size` of the kept cells, as the Signac vignettes do, or 1e4
  on an ATAC-only experiment, whose `lib_size` counts peaks.

- barcode_map:

  Optional. Cell id to fragment barcode, where they differ: a named
  character vector or a function.

- overwrite:

  Boolean. Replace an existing store at `dir_out`.

- max_mem_gb:

  Optional numeric. Memory budget of the gene-major conversion.

- .verbose:

  Boolean. Controls verbosity.

## Value

The gene activity experiment as a loaded `SingleCells`. var holds
`gene_id` (the GTF gene ID without its version, as Cell Ranger keys the
RNA), `gene_name` (made unique), `chrom`, `start` and `end`; obs's
`lib_size` and `nnz` are the gene activity's.

## References

Stuart, et al., Nat Methods, 2021
