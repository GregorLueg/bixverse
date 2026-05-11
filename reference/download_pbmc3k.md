# Download PBMC3K data from Zenodo

This function downloads the PBMC3K dataset from 10x Genomics (uploaded
on Zenodo) and extracts it and returns the (temporary) paths.

## Usage

``` r
download_pbmc3k(quiet = FALSE)
```

## Arguments

- quiet:

  Boolean. If the download shall be quiet.

## Value

String. The path to the extracted PBMC3K data.
