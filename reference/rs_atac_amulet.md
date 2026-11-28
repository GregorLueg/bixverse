# AMULET multiplet detection from fragment overlaps

**\[experimental\]**

## Usage

``` r
rs_atac_amulet(
  fragments,
  barcodes,
  restrict_beds,
  max_overlap,
  max_fragment_len,
  repeat_q,
  multiplet_q,
  chroms,
  max_mem_gb,
  count_trailing_loci
)
```

## Arguments

- fragments:

  String. The fragment file.

- barcodes:

  Character vector. Cells to test.

- restrict_beds:

  Character vector. BED files of regions to drop (repeats, segmental
  duplications, blacklists).

- max_overlap:

  Integer. Overlapping fragments a diploid cell can explain.

- max_fragment_len:

  Integer. Longest valid fragment.

- repeat_q:

  Numeric. q-value under which a site is a repeat.

- multiplet_q:

  Numeric. q-value under which a cell is a multiplet.

- chroms:

  Optional character vector. Chromosomes to use, `NULL` for the
  autosomes.

- max_mem_gb:

  Optional numeric. Read budget, `NULL` for one chromosome at a time.

- count_trailing_loci:

  Boolean. Count loci no later fragment passes.

## Value

A list with per-cell `n_overlaps`, `n_sites`, `p_value`, `q_value` and
`multiplet`.

## References

Thibodeau, et al., Genome Biol, 2021
