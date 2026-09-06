# Read in the motif rankings and transform them into a matrix

This function loads in the .feather files with the motif to target gene
rankings. These can be found here:
`https://resources.aertslab.org/cistarget/databases/`

## Usage

``` r
read_motif_ranking(ranking_file)
```

## Arguments

- ranking_file:

  String. The file path to the .feather file

## Value

An integer matrix that has been transposed for easier use in the
underlying Rust code.

## Examples

``` r
if (FALSE) { # \dontrun{
# transposed motif rankings from the downloaded feather file
paths <- download_cistarget_hg38()
rankings <- read_motif_ranking(paths$rankings)
dim(rankings)
} # }
```
