# Count a fragment file into peaks or genome tiles as the ATAC store

Counts every cell of the RNA store, in its order, so the RNA cell filter
applies unchanged. Pass exactly one of `peaks` and `tile_width`. The
features go to `var_atac` with `chr:start-end` ids (0-based, half-open),
the per-cell totals to obs as `atac_lib_size` and `atac_nnz`.

## Usage

``` r
add_atac_fragments_sc(
  object,
  fragments,
  peaks = NULL,
  tile_width = NULL,
  chrom_sizes = NULL,
  mode = c("insertion", "paired-insertion", "fragment"),
  barcode_map = NULL,
  overwrite = FALSE,
  max_mem_gb = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCellsMultiModal` class.

- fragments:

  String. The fragment file (BGZF, gzip or plain). With a tabix index
  next to it, passes read one chromosome at a time.

- peaks:

  Optional. A BED file or a data.frame with `chrom`, `start`, `end`
  (0-based, half-open), e.g. from
  [`call_peaks_sc()`](https://gregorlueg.github.io/bixverse/reference/call_peaks_sc.md).

- tile_width:

  Optional integer. Tile width in bp, instead of `peaks`.

- chrom_sizes:

  Optional named integer vector of chromosome lengths, for the tiles.

- mode:

  String. What one count is. `"insertion"` (each Tn5 insertion, Cell
  Ranger and ArchR), `"paired-insertion"` (once per fragment and
  feature, SnapATAC2) or `"fragment"` (each overlapping fragment,
  Signac).

- barcode_map:

  Optional. Cell id to fragment barcode, where they differ: a named
  character vector (names are the cell ids) or a function.

- overwrite:

  Boolean. Replace an existing ATAC store.

- max_mem_gb:

  Optional numeric. Memory budget of the peak-major conversion.

- .verbose:

  Boolean. Controls verbosity.

## Value

The object with the ATAC store attached.
