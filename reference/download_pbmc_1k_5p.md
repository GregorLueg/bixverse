# Download the raw PBMC 1k 5' matrix from 10x Genomics

Downloads the unfiltered `raw_feature_bc_matrix.h5` of the 10x Genomics
5' PBMC 1k run, straight from the 10x CDN rather than the bixverse-data
mirror. All 737,280 barcodes are in there, empty droplets included,
which is what
[`cellsweep_sc()`](https://gregorlueg.github.io/bixverse/reference/cellsweep_sc.md)
needs. The file also carries 19 Antibody Capture features next to the
36,601 genes, so load it with
`load_tenx_h5(feature_type = "Gene Expression")`.

## Usage

``` r
download_pbmc_1k_5p(quiet = FALSE)
```

## Arguments

- quiet:

  Boolean. If the download shall be quiet.

## Value

String. The path to the downloaded h5 file.

## Examples

``` r
if (FALSE) { # \dontrun{
# pulls roughly 17MB into the session tempdir()
h5_path <- download_pbmc_1k_5p()
read_tenx_h5_metadata(h5_path)$dims
} # }
```
