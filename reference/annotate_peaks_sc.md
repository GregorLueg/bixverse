# The closest gene of every ATAC feature

Signac's `ClosestFeature()`. Written to the ATAC feature table as
`closest_gene`, `closest_gene_id`, `distance_to_gene` (bases to the gene
body, 0 inside or touching it) and `distance_to_tss` (feature centre to
TSS, positive downstream).

## Usage

``` r
annotate_peaks_sc(object, gtf, biotypes = "protein_coding")
```

## Arguments

- object:

  `SingleCellsAtac` or `SingleCellsMultiModal` class.

- gtf:

  String. Gene annotation.

- biotypes:

  Optional character vector. Biotypes to consider. `NULL` keeps all.

## Value

The object with the four columns added to the ATAC feature table.
