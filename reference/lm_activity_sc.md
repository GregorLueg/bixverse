# Calculate ULM / MLM activity scores per cell

Rust implementation of decoupler's univariate (ULM) and multivariate
(MLM) linear model activity scores, per cell (or meta cell) on the
normalised counts. See
[`calc_lm_activity()`](https://gregorlueg.github.io/bixverse/reference/calc_lm_activity.md)
for the model. Cells are read and scored in chunks, so memory stays
bounded on large data.

Differences to decoupler and things to know:

- Every gene of the store enters the model, including genes that are
  zero in every cell. decoupler drops those, so scores differ slightly
  where targets were never detected.

- Empty cells score 0 with a p-value of 1.

- MLM errors on collinear sources.

- The p-values treat genes as independent, which they are not in a
  single cell. Read them as descriptive.

## Usage

``` r
lm_activity_sc(object, net, lm_params = params_lm_activity(), .verbose = TRUE)
```

## Arguments

- object:

  `SingleCells` or `MetaCells` class.

- net:

  data.frame. The network in long format with the columns `source`,
  `target` and optionally `weight`. Targets need to be gene identifiers
  of the object.

- lm_params:

  List. The parameters, see
  [`params_lm_activity()`](https://gregorlueg.github.io/bixverse/reference/params_lm_activity.md).

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

For `SingleCells`, an `ScMatrixRes` of cells x sources with the
activities if `p_values = FALSE`. Otherwise a list with the
`ScMatrixRes` elements `es`, `pvals` and `fdr` (BH per cell across
sources). For `MetaCells` the same, but plain matrices of meta cells x
sources. `NULL` with a warning if no source has at least `tmin` targets
in the object.

## References

Badia-i-Mompel, et al., Bioinform. Adv., 2022

## Examples

``` r
# score two made-up regulons on the demo data
sc <- demo_single_cells()
genes <- get_gene_names(sc)
net <- data.frame(
  source = rep(c("tf_a", "tf_b"), each = 10L),
  target = genes[1:20],
  weight = 1
)
res <- lm_activity_sc(sc, net = net, .verbose = FALSE)
dim(res)
#> [1] 500   2

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
