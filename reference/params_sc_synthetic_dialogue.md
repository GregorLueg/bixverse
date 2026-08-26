# Default parameters for generation of synthetic DIALOGUE data

Shapes the fixture with a planted multicellular programme that
[`generate_dialogue_test_data()`](https://gregorlueg.github.io/bixverse/reference/generate_dialogue_test_data.md)
builds. The defaults give 14 samples x 25 cells x 3 cell types = 1050
cells over 400 genes.

## Usage

``` r
params_sc_synthetic_dialogue(
  n_samples = 14L,
  cells_per_sample = 25L,
  n_cell_types = 3L,
  n_features = 8L,
  n_sample_features = 5L,
  n_genes = 400L,
  n_planted = 8L
)
```

## Arguments

- n_samples:

  Integer. Samples the experiment spans. DIALOGUE needs at least 5.

- cells_per_sample:

  Integer. Cells per sample per cell type.

- n_cell_types:

  Integer. Number of cell types. Must be at least 2.

- n_features:

  Integer. Feature columns per cell type. Must be at least 2.

- n_sample_features:

  Integer. Feature columns carrying a per-sample component. The first of
  those is the shared programme, the rest are cell-type-specific
  nuisance; anything past this count is pure noise and exists so the
  ANOVA filter has something to reject.

- n_genes:

  Integer. Number of genes.

- n_planted:

  Integer. Planted genes per cell type. The blocks are contiguous, so
  `n_planted * n_cell_types` has to fit into `n_genes`.

## Value

A list with the parameters.

## Details

`n_genes` defaults higher here than on the Rust side, which plants the
same structure into 90. R runs the counts through the normal ingestion
path, and on a small panel the planted genes are a large enough share of
the library that log-normalisation divides the signal back out: at 90
genes the library size tracks the programme at an r of 0.75 and
background genes pick up a spurious correlation of their own. At 400 the
planted block is a few percent of the library and the contrast survives.
The Rust tests feed the planted layer straight in, so they never meet
this.

## References

Jerby-Arnon & Regev, Nature Biotechnology, 2022
