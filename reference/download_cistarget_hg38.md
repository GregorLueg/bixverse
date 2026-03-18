# Download CisTarget reference files for human (hg38)

Download CisTarget reference files for human (hg38)

## Usage

``` r
download_cistarget_hg38(
  cache_dir = tools::R_user_dir("bixverse", which = "cache"),
  overwrite = FALSE
)
```

## Arguments

- cache_dir:

  String. Directory to cache the files. Defaults to a package-specific
  user cache directory.

- overwrite:

  Logical. Re-download even if files already exist.

## Value

Named list with paths: `rankings` and `motif_annotations`.
