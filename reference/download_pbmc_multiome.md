# Download the 10x PBMC unsorted 3k multiome data

The 10x PBMC unsorted 3k multiome (Cell Ranger ARC 2.0): the filtered
feature barcode matrix (RNA and Cell Ranger's peak counts in one h5),
the fragment file with its tabix index, Cell Ranger's peaks, the GENCODE
v32 annotation its GRCh38 reference was built from, and the ENCODE hg38
blacklist (v2, Amemiya et al. 2019). 10x Genomics data, CC BY 4.0. Files
already in `dir_out` are not downloaded again.

## Usage

``` r
download_pbmc_multiome(dir_out = NULL, quiet = FALSE)
```

## Arguments

- dir_out:

  Optional string. Where to put the files. `NULL` uses `pbmc_multiome`
  in the session [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- quiet:

  Boolean. If the download shall be quiet.

## Value

String. The directory with the files under their 10x names.

## Examples

``` r
if (FALSE) { # \dontrun{
# pulls the files into the session tempdir()
path <- download_pbmc_multiome()
list.files(path)
} # }
```
