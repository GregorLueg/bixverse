# Download the marrow CD34 example data from Palantir

This function downloads the bone marrow CD34 data set from the Palantir
paper into the temporary directory.

## Usage

``` r
download_marrow_cd34(quiet = FALSE)
```

## Arguments

- quiet:

  Boolean. If the download shall be quiet.

## Value

String. The path to the marrow CD34 data set.

## References

Setty, et al., Nat. Biotechnol., 2019

## Examples

``` r
if (FALSE) { # \dontrun{
# pulls the archive into the session tempdir()
path <- download_marrow_cd34()
get_h5ad_dimensions(path)$dims
} # }
```
