# Resolve the download URLs for an example data file

GitHub releases first, Zenodo second. Zenodo holds the archival copy
with the DOI, but it is slow and drops connections often enough to fail
a vignette build outright.

## Usage

``` r
.data_urls(file, zenodo_record = NULL)
```

## Arguments

- file:

  String. File name, identical in both locations.

- zenodo_record:

  Optional string. Zenodo record holding the fallback copy. `NULL` for
  datasets that only live in the GitHub release, which is the case until
  they are archived.

## Value

Character vector of URLs, in the order they should be tried.
