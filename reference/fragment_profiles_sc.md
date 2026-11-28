# TSS insertion profile and fragment length histogram of the kept cells

One pass over the fragment file, per level of `group_by` (SnapATAC2's
`TSS_profile` and `frag_size_distr`, which these reproduce): insertions
are a fragment's start and `end - 1`, offsets are strand aware, positive
downstream.

## Usage

``` r
fragment_profiles_sc(
  object,
  fragments,
  gtf,
  group_by = NULL,
  half_width = 2000L,
  max_len = 1000L,
  flank = 100L,
  biotypes = "protein_coding",
  barcode_map = NULL
)
```

## Arguments

- object:

  `SingleCellsAtac` or `SingleCellsMultiModal` class.

- fragments:

  String. The fragment file.

- gtf:

  String. Gene annotation for the TSSs.

- group_by:

  Optional string. obs column to split by; `NULL` for one group,
  `"all"`.

- half_width:

  Integer. Profile half width in bp.

- max_len:

  Integer. Longest fragment length in the histogram.

- flank:

  Integer. Bases at each end of the profile whose mean is the background
  of the enrichment.

- biotypes:

  Optional character vector. Gene biotypes whose TSSs count. `NULL`
  keeps all.

- barcode_map:

  Optional. Cell id to fragment barcode, where they differ.

## Value

A list with

- tss - data.table with `group`, `offset`, `insertions` and `enrichment`
  (insertions over the mean of both flanks).

- lengths - data.table with `group`, `length`, `fragments` and
  `fraction` (of the group's fragments up to `max_len`).
