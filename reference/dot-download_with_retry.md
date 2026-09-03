# Download a file, retrying on failure

Zenodo drops connections often enough that a single
[`download.file()`](https://rdrr.io/r/utils/download.file.html) call is
not reliable, which shows up as `cannot open URL` mid-vignette. Retries
with a linear backoff.

## Usage

``` r
.download_with_retry(urls, dest_file, quiet = FALSE, tries = 4L)
```

## Arguments

- urls:

  Character vector. Mirrors to try, in order.

- dest_file:

  String. Path to write the file to.

- quiet:

  Boolean. If the download shall be quiet.

- tries:

  Integer. Attempts before giving up.

## Value

Invisibly, `dest_file`.
