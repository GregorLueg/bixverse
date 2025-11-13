# Read in the motif annotation file

This function loads in the motif2tf information that you can get from
`https://resources.aertslab.org/cistarget/motif2tf/`. The function will
generate a data.table that can be subsequently used.

## Usage

``` r
read_motif_annotation_file(annot_file)
```

## Arguments

- annot_file:

  String. Path to the motif2tf file that you downloaded.

## Value

data.table with the motif to transcription factor information.
