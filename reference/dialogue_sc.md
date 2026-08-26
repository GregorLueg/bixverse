# Find multicellular programmes with DIALOGUE

DIALOGUE looks for programmes of cell-type-specific genes whose activity
covaries across the samples that several cell types share. Where
co-expression modules ask what varies together *within* a cell type,
this asks what varies together *between* them, sample by sample.

## Usage

``` r
dialogue_sc(
  object,
  cell_type_col,
  sample_col,
  features,
  quality_col = NULL,
  gene_ids = NULL,
  pmd_params = params_dialogue_pmd(),
  hlm_params = params_dialogue_hlm(),
  refine_params = params_dialogue_refine(),
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells`, `SingleCellsSubset` or `MetaCells` class.

- cell_type_col:

  String. Column in the obs table holding the cell type labels.

- sample_col:

  String. Column in the obs table holding the sample labels. The random
  intercept in stage two is over these. For `MetaCells` the meta cells
  must have been built *within* samples, otherwise the level is not
  well-defined.

- features:

  Named list of numeric matrices, one per cell type. Names must match
  the levels in `cell_type_col`, row names must cover that cell type's
  cells, and each needs at least two columns. Rows are matched by name,
  not by position.

- quality_col:

  Optional string. Column in the obs table to use as the cell quality
  covariate. If `NULL`, defaults to the z-scored log library size.

- gene_ids:

  Optional character. Genes to consider when building signatures. If
  `NULL`, uses
  [`get_hvg()`](https://gregorlueg.github.io/bixverse/reference/get_hvg.md)
  on the object.

- pmd_params:

  List, see
  [`params_dialogue_pmd()`](https://gregorlueg.github.io/bixverse/reference/params_dialogue_pmd.md).

- hlm_params:

  List, see
  [`params_dialogue_hlm()`](https://gregorlueg.github.io/bixverse/reference/params_dialogue_hlm.md).

- refine_params:

  List, see
  [`params_dialogue_refine()`](https://gregorlueg.github.io/bixverse/reference/params_dialogue_refine.md).

- .verbose:

  Boolean or integer. Verbosity.

## Value

A `DialogueResult` object.

## Details

The algorithm works in three stages. Each cell type's features are
collapsed to one row per sample and put through a sparse multi-CCA,
giving every programme a weight vector per cell type plus a provisional
gene signature. Then, for every ordered pair of cell types and every
candidate gene, a mixed model asks whether a cell's own programme score
tracks the partner's expression of that gene in the same sample. Finally
the partners are meta-analysed and the scores refit onto the surviving
genes by non-negative least squares.

`features` is mandatory and there is no default. DIALOGUE does not
compute it: it is whatever low-dimensional description of each cell type
you trust, and the method is only as good as that choice. Do *not* hand
it a slice of a global PCA. Those components mostly carry
between-cell-type identity, which is near-constant inside one cell type
and gets dropped by the ANOVA filter, leaving the decomposition to work
off whatever is left. Run a PCA per cell type instead, which for
`SingleCells` means
[`SingleCellsSubset()`](https://gregorlueg.github.io/bixverse/reference/SingleCellsSubset.md)
followed by
[`calculate_pca_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_pca_sc.md)
on each subset.

The method is unforgiving about study design, and the failure modes are
errors rather than bad answers. It needs at least two cell types, at
least five samples present in *every* cell type, and enough cells per
sample per cell type to clear `abn_c` in
[`params_dialogue_pmd()`](https://gregorlueg.github.io/bixverse/reference/params_dialogue_pmd.md).
A handful of samples with thousands of cells each is the regime it was
built for; many samples with a dozen cells each is not.

## References

Jerby-Arnon & Regev, Nature Biotechnology, 2022
