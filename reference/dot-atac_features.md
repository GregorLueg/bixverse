# The features to count fragments into

The features to count fragments into

## Usage

``` r
.atac_features(peaks, tile_width, chrom_sizes)
```

## Arguments

- peaks:

  `NULL`, a BED file or a data.frame of regions.

- tile_width:

  Optional integer. Tile width in bp, instead of `peaks`.

- chrom_sizes:

  Optional named integer vector of chromosome lengths, for the tiles.

## Value

A list with `chrom`, `start`, `end` and `kind` (`"peak"` or `"tile"`).
