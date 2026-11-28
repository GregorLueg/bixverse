# Load an ATAC-only experiment from a fragment file

The cells come from exactly one of `barcodes`, `cellranger` or
`min_fragments`; with the last, every barcode with at least that many
fragments (outside `exclude_chroms`), sorted, as Cell Ranger orders
them. See
[`barcode_fragment_counts()`](https://gregorlueg.github.io/bixverse/reference/barcode_fragment_counts.md)
for the knee plot to pick it from. Pass exactly one of `peaks` and
`tile_width`.

## Usage

``` r
load_atac_fragments(
  object,
  fragments,
  peaks = NULL,
  tile_width = NULL,
  chrom_sizes = NULL,
  mode = c("insertion", "paired-insertion", "fragment"),
  barcodes = NULL,
  cellranger = NULL,
  min_fragments = NULL,
  exclude_chroms = c("chrM", "MT", "chrMT", "M"),
  max_mem_gb = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCellsAtac` class.

- fragments:

  String. The fragment file (BGZF, gzip or plain).

- peaks:

  Optional. A BED file or a data.frame with `chrom`, `start`, `end`.

- tile_width:

  Optional integer. Tile width in bp, instead of `peaks`.

- chrom_sizes:

  Optional named integer vector of chromosome lengths, for the tiles.

- mode:

  String. What one count is, see
  [`add_atac_fragments_sc()`](https://gregorlueg.github.io/bixverse/reference/add_atac_fragments_sc.md).

- barcodes:

  Optional character vector. The cell barcodes.

- cellranger:

  Optional string. Cell Ranger ATAC's `singlecell.csv` or filtered h5.

- min_fragments:

  Optional integer. Call cells with at least this many fragments.

- exclude_chroms:

  Character vector. Chromosomes the `min_fragments` count skips.

- max_mem_gb:

  Optional numeric. Memory budget of the feature-major conversion.

- .verbose:

  Boolean. Controls verbosity.

## Value

The loaded experiment. obs holds `cell_id`, `lib_size`, `nnz`,
`to_keep`; var the features as `gene_id` (`chr:start-end`) with `chrom`,
`start`, `end` and `no_cells_exp`.
