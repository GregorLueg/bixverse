# Flag multiplets from fragment overlaps (AMULET)

A diploid nucleus yields at most two overlapping fragments at a locus.
Per cell, the loci with more are counted (autosomes, restriction regions
dropped), sites shared by too many cells are dropped as repeats, and the
count is tested against the mean (Poisson, Benjamini-Hochberg).
Implemented from the paper; on the test data it reproduces the reference
implementation's counts and p-values. Needs no clustering and no RNA, so
it works for ATAC-only experiments. Run it per library: multiplets are a
property of one capture.

Adds `amulet_n_overlaps`, `amulet_n_sites`, `amulet_p`, `amulet_q` and
`amulet_multiplet` to obs for every cell of the store.

## Usage

``` r
amulet_sc(
  object,
  fragments,
  restrict_beds = character(),
  amulet_params = params_atac_amulet(),
  barcode_map = NULL,
  max_mem_gb = NULL
)
```

## Arguments

- object:

  `SingleCells` class; every cell of its store is tested.

- fragments:

  String. The fragment file.

- restrict_beds:

  Character vector. BED files of regions to ignore: repeats, segmental
  duplications, blacklists. The paper uses all three.

- amulet_params:

  List. Output of
  [`params_atac_amulet()`](https://gregorlueg.github.io/bixverse/reference/params_atac_amulet.md).

- barcode_map:

  Optional. Cell id to fragment barcode, where they differ.

- max_mem_gb:

  Optional numeric. Read budget; `NULL` reads one chromosome at a time
  (with a tabix index).

## Value

The object with the AMULET columns added to obs.

## References

Thibodeau, et al., Genome Biol, 2021
