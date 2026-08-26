# Run DIALOGUE over a set of single cells

**\[experimental\]** Finds multicellular programmes: genes that are
cell-type-specific but whose activity covaries across the samples the
cell types share. Three stages. First, each cell type's features are
collapsed to one row per sample and put through a sparse multi-CCA,
giving every programme a weight vector per cell type plus a provisional
gene signature. Second, for every ordered pair of cell types and every
candidate gene, a mixed model asks whether a cell's own programme score
tracks the partner cell type's expression of that gene in the same
sample. Third, the partners are meta-analysed and the scores refit onto
the surviving genes by non-negative least squares.

The counts are streamed off the gene-major file and only the normalised
layer is read. DIALOGUE does not compute the features: whatever the
caller trusts as a low-dimensional description of each cell type goes in
as `features`.

## Usage

``` r
rs_dialogue_sc(
  f_path_gene,
  cell_type_indices,
  features,
  sample_ids,
  cell_quality,
  gene_indices,
  dialogue_params,
  verbose
)
```

## Arguments

- f_path_gene:

  Path to the `counts_genes.bin` file.

- cell_type_indices:

  List of integer vectors. 0-indexed(!) positions of the cells belonging
  to each cell type. At least two cell types are needed.

- features:

  List of numeric matrices, one per cell type, shaped
  `n_cells_in_type x k_i` with rows aligned to `cell_type_indices`.
  Needs at least two columns per cell type.

- sample_ids:

  Integer vector. 0-indexed(!) sample code per cell, over the *whole*
  store rather than per cell type. Must be long enough to cover the
  largest index in `cell_type_indices`.

- cell_quality:

  Numeric vector. Quality covariate per cell, indexed the same way as
  `sample_ids`. Upstream's `cellQ`, typically the z-scored log of the
  library size.

- gene_indices:

  Integer vector. 0-indexed(!) positions of the genes to consider when
  building signatures.

- dialogue_params:

  Named list. Contains the DIALOGUE parameters across all three stages,
  see
  [`params_dialogue_pmd()`](https://gregorlueg.github.io/bixverse/reference/params_dialogue_pmd.md),
  [`params_dialogue_hlm()`](https://gregorlueg.github.io/bixverse/reference/params_dialogue_hlm.md)
  and
  [`params_dialogue_refine()`](https://gregorlueg.github.io/bixverse/reference/params_dialogue_refine.md).
  The three blocks share one flat list.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- shared_samples - Integer vector. 0-indexed(!) sample codes present in
  every cell type.

- kept_features - List of integer vectors. 0-indexed(!) feature columns
  surviving the ANOVA filter, per cell type.

- mcp_cell_types - List of integer vectors. 0-indexed(!) cell types each
  programme spans.

- ws - List of matrices. Sparse canonical weights per cell type,
  `kept_features x k`.

- scores - List of matrices. Final programme scores per cell type,
  `n_cells_in_type x k`, residualised on the quality covariate.

- cca_scores - List of matrices. Stage one's canonical scores, kept for
  comparison against the refit.

- emp_p - Matrix. Empirical p per programme and cell type pair,
  `k x n_pairs`, the columns in `combn(n_cell_types, 2)` order.

- pair_cor - Matrix. Canonical correlation on the real fit, same shape.

- refit_fidelity - Matrix. Correlation between the canonical score and
  the refit, `n_cell_types x k`. A low value means the gene-level refit
  drifted away from the programme the decomposition found.

- verdicts - List of equal-length vectors with the meta-analysis verdict
  per gene: `cell_type`, `programme`, `gene` (all 0-indexed), `up`,
  `n_supporting`, `support_fraction`, `p_up`, `p_down` and
  `coefficient`.

- permissive - Nested list `[cell_type][programme]` of `up` / `down`
  gene positions (0-indexed).

- strict - The same, for the stricter gene list.

## References

Jerby-Arnon & Regev, Nature Biotechnology, 2022
