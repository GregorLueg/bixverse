# Convert legacy v2 single-cell data files to v3 format

Takes a directory containing legacy `counts_cells.bin` and
`counts_genes.bin` files and converts them to the v3 format in place.

## Usage

``` r
sc_old_file_conversion(dir_path, verbose = TRUE)
```

## Arguments

- dir_path:

  Character scalar. Path to the directory containing the legacy binary
  files.

- verbose:

  Logical scalar. Whether to print progress messages. Default `TRUE`.

## Value

Invisible `NULL`. Called for its side effects.
