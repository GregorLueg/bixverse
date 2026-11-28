# Check a fragment file and note what it lacks for a fast read

Check a fragment file and note what it lacks for a fast read

## Usage

``` r
.check_fragments(fragments, .verbose = TRUE)
```

## Arguments

- fragments:

  String. The fragment file.

- .verbose:

  Boolean. Controls verbosity.

## Value

Invisible `NULL`. Errors when the file is missing, warns when the tabix
index is older than the file.
