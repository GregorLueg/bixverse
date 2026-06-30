# Bixverse implementation of singscore (multiple gene sets)

Scores all samples against many up-regulated gene sets with optional
paired down-regulated sets. Down sets are paired with up sets by name;
unmatched names are dropped.

## Usage

``` r
calc_singscore_multi(
  ranks,
  up_pathways,
  down_pathways = NULL,
  center_score = TRUE,
  known_direction = TRUE,
  min_size = 1L,
  max_size = 500L
)
```

## Arguments

- ranks:

  Numerical matrix. Output of
  [`calc_singscore_rank()`](https://gregorlueg.github.io/bixverse/reference/calc_singscore_rank.md).

- up_pathways:

  Named list of character vectors. Up-regulated gene sets.

- down_pathways:

  Named list or NULL. Paired down-regulated gene sets.

- center_score:

  Boolean. Centre scores around 0. Ignored when
  `known_direction = FALSE`.

- known_direction:

  Boolean. Whether the up-set direction is known. Becomes irrelevant
  when `down_set` is also provided.

- min_size:

  Integer. Minimum gene-set size after dropping missing genes.

- max_size:

  Integer. Maximum gene-set size.

## Value

A named list with two matrices, `scores` and `dispersions`, each of
shape gene_sets × samples.

## References

Foroutan et al., BMC Bioinformatics, 2018.
