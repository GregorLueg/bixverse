# Bixverse implementation of singscore (single gene set)

Scores all samples against one up-regulated gene set with an optional
paired down-regulated set. When `n_permutations > 0`, additionally runs
a permutation test and includes empirical p-values.

## Usage

``` r
calc_singscore(
  ranks,
  up_set,
  down_set = NULL,
  center_score = TRUE,
  known_direction = TRUE,
  n_permutations = 0L,
  seed = 42L
)
```

## Arguments

- ranks:

  Numerical matrix. Output of
  [`calc_singscore_rank()`](https://gregorlueg.github.io/bixverse/reference/calc_singscore_rank.md).

- up_set:

  Character vector. Gene names of the up-regulated set.

- down_set:

  Character vector or NULL. Optional paired down-regulated set.

- center_score:

  Boolean. Centre scores around 0. Ignored when
  `known_direction = FALSE`.

- known_direction:

  Boolean. Whether the up-set direction is known. Becomes irrelevant
  when `down_set` is also provided.

- n_permutations:

  Integer. Number of permutations. `0` disables the permutation test.

- seed:

  Integer. RNG seed for the permutation test.

## Value

A data.table with one row per sample. Columns: `total_score`,
`total_dispersion`, optionally `up_score`, `up_dispersion`,
`down_score`, `down_dispersion`, and (when permutations are run) `pval`.
With permutations, the null distribution is attached as
`attr(., "null_distribution")`.

## References

Foroutan et al., BMC Bioinformatics, 2018.

## Examples

``` r
# score every sample against one up-regulated gene set
set.seed(123L)
exp_mat <- matrix(
  rnorm(200 * 10),
  nrow = 200,
  dimnames = list(sprintf("gene_%03i", 1:200), sprintf("sample_%i", 1:10))
)
ranks <- calc_singscore_rank(exp_mat)
head(calc_singscore(ranks, up_set = sprintf("gene_%03i", 1:20)))
#>    total_score total_dispersion sample_id
#>          <num>            <num>    <char>
#> 1:  0.05611111         78.57792  sample_1
#> 2: -0.01944444         61.52799  sample_2
#> 3: -0.07138889         57.08019  sample_3
#> 4: -0.08083333         49.66717  sample_4
#> 5:  0.01888889         57.82149  sample_5
#> 6:  0.03388889         74.87141  sample_6
```
