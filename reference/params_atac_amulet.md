# Wrapper function for the AMULET parameters

Parameters of
[`amulet_sc()`](https://gregorlueg.github.io/bixverse/reference/amulet_sc.md),
multiplet detection from fragment overlaps.

## Usage

``` r
params_atac_amulet(
  max_overlap = 2L,
  max_fragment_len = 900L,
  repeat_q = 0.01,
  multiplet_q = 0.01,
  chroms = NULL,
  count_trailing_loci = FALSE
)
```

## Arguments

- max_overlap:

  Integer. Overlapping fragments a cell can explain; 2 for diploid.
  Defaults to `2L`.

- max_fragment_len:

  Integer. Longest valid fragment. Defaults to `900L`.

- repeat_q:

  Numeric. q-value under which a site is a repeat. Defaults to `0.01`.

- multiplet_q:

  Numeric. q-value under which a cell is a multiplet. Defaults to
  `0.01`.

- chroms:

  String vector or `NULL`. Chromosomes to use. `NULL` uses the
  autosomes. Defaults to `NULL`.

- count_trailing_loci:

  Boolean. Count a cell's last loci on the last chromosome used, which
  the reference implementation drops. Defaults to `FALSE`.

## Value

A named list with the following elements:

- max_overlap - Integer. Overlapping fragments a cell can explain; 2 for
  diploid. Defaults to `2L`.

- max_fragment_len - Integer. Longest valid fragment. Defaults to
  `900L`.

- repeat_q - Numeric. q-value under which a site is a repeat. Defaults
  to `0.01`.

- multiplet_q - Numeric. q-value under which a cell is a multiplet.
  Defaults to `0.01`.

- chroms - String vector or `NULL`. Chromosomes to use. `NULL` uses the
  autosomes. Defaults to `NULL`.

- count_trailing_loci - Boolean. Count a cell's last loci on the last
  chromosome used, which the reference implementation drops. Defaults to
  `FALSE`.

## References

Thibodeau, et al., Genome Biol, 2021.
