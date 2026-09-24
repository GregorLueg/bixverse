# Wrapper function for the DIALOGUE decomposition parameters

Stage one of DIALOGUE: the penalised matrix decomposition that turns the
per-cell-type features into multicellular programmes, and the
provisional gene signatures that come off it.

## Usage

``` r
params_dialogue_pmd(
  k = 2L,
  n_permutations = 100L,
  extra_sparse = FALSE,
  abn_c = 15L,
  p_anova = 0.05,
  centre = TRUE,
  cap = 0.01,
  spatial = FALSE,
  n_genes = 200L,
  min_ci = 0.05,
  averaging = c("median", "mean"),
  mcp_assignment_p = 0.1,
  seed = 1234L
)
```

## Arguments

- k:

  Integer. Number of multicellular programmes to extract. Must be at
  least 1. Defaults to `2L`.

- n_permutations:

  Integer. Permutations backing the empirical p-value per programme.
  Must be at least 2. Defaults to `100L`.

- extra_sparse:

  Boolean. Tune the L1 bound by permutation instead of fixing it at
  `sqrt(p_1) / 2`. Costs ten more fits per permutation. Defaults to
  `FALSE`.

- abn_c:

  Integer. Minimum cells a sample must contribute, within a cell type,
  before it counts towards the feature-level ANOVA. Defaults to `15L`.

- p_anova:

  Numeric. BH-adjusted ANOVA cutoff for keeping a feature. Must be in
  `(0, 1]`. Defaults to `0.05`.

- centre:

  Boolean. Centre and scale the sample-level feature matrix, then
  winsorise it. Defaults to `TRUE`.

- cap:

  Numeric. Winsorising tail fraction applied to each column. Must be in
  `[0, 0.5)`. Defaults to `0.01`.

- spatial:

  Boolean. Spatial data: skip the ANOVA feature filter entirely. Niches
  are small, so a feature need not vary across them to be real. Defaults
  to `FALSE`.

- n_genes:

  Integer. Genes taken per programme per direction when building a
  signature. Defaults to `200L`.

- min_ci:

  Numeric. Minimum absolute correlation for a gene to enter a signature.
  Must be in `[0, 1]`. Defaults to `0.05`.

- averaging:

  String. How cell-level features are collapsed per sample. One of
  `c("median", "mean")`. Defaults to `"median"`.

- mcp_assignment_p:

  Numeric. Empirical p below which a cell type pair counts as connected
  when deciding which cell types a programme spans. Must be in `(0, 1]`.
  Defaults to `0.1`.

- seed:

  Integer. Seed for the permutation null. Defaults to `1234L`.

## Value

A named list with the following elements:

- k - Integer. Number of multicellular programmes to extract. Must be at
  least 1. Defaults to `2L`.

- n_permutations - Integer. Permutations backing the empirical p-value
  per programme. Must be at least 2. Defaults to `100L`.

- extra_sparse - Boolean. Tune the L1 bound by permutation instead of
  fixing it at `sqrt(p_1) / 2`. Costs ten more fits per permutation.
  Defaults to `FALSE`.

- abn_c - Integer. Minimum cells a sample must contribute, within a cell
  type, before it counts towards the feature-level ANOVA. Defaults to
  `15L`.

- p_anova - Numeric. BH-adjusted ANOVA cutoff for keeping a feature.
  Must be in `(0, 1]`. Defaults to `0.05`.

- centre - Boolean. Centre and scale the sample-level feature matrix,
  then winsorise it. Defaults to `TRUE`.

- cap - Numeric. Winsorising tail fraction applied to each column. Must
  be in `[0, 0.5)`. Defaults to `0.01`.

- spatial - Boolean. Spatial data: skip the ANOVA feature filter
  entirely. Niches are small, so a feature need not vary across them to
  be real. Defaults to `FALSE`.

- n_genes - Integer. Genes taken per programme per direction when
  building a signature. Defaults to `200L`.

- min_ci - Numeric. Minimum absolute correlation for a gene to enter a
  signature. Must be in `[0, 1]`. Defaults to `0.05`.

- averaging - String. How cell-level features are collapsed per sample.
  One of `c("median", "mean")`. Defaults to `"median"`.

- mcp_assignment_p - Numeric. Empirical p below which a cell type pair
  counts as connected when deciding which cell types a programme spans.
  Must be in `(0, 1]`. Defaults to `0.1`.

- seed - Integer. Seed for the permutation null. Defaults to `1234L`.

## Details

The defaults follow upstream's `DLG.get.param`. Two knobs are worth
thinking about before anything else. `k` is how many programmes you are
asking for, and there is no sweep to help you pick it. `n_permutations`
sets the resolution of the empirical p-value: with the default of `100`
the smallest p you can observe is `0.01`, so lower it for a quick look
and leave it alone for anything you intend to believe.

`averaging` is exposed and honoured here. Upstream takes the same
argument and then ignores it, hard-coding column medians, so `"median"`
is what every published DIALOGUE run actually used.

## References

Jerby-Arnon & Regev, Nature Biotechnology, 2022
