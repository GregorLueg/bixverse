# Wrapper function for the DIALOGUE mixed model parameters

Stage two of DIALOGUE: for every ordered pair of cell types and every
candidate gene, a random-intercept mixed model over samples asking
whether a cell's own programme score tracks the partner cell type's
expression of that gene in the same sample.

## Usage

``` r
params_dialogue_hlm(
  min_cells_per_sample = 2L,
  use_tme_qc = TRUE,
  use_cell_quality = TRUE,
  satterthwaite = TRUE
)
```

## Arguments

- min_cells_per_sample:

  Integer. Minimum cells a sample must contribute, in *both* cell types
  of a pair, before it takes part in that pair's models.

- use_tme_qc:

  Boolean. Include the partner cell type's mean quality in that sample
  as a fixed effect. Upstream's `tme.qc`.

- use_cell_quality:

  Boolean. Include the responding cell's own quality as a fixed effect.
  Upstream's `cellQ`.

- satterthwaite:

  Boolean. Compute Satterthwaite denominator degrees of freedom, as
  `lmerTest` does.

## Value

A list with the stage two DIALOGUE parameters.

## Details

This stage dominates the runtime, and `satterthwaite` is the knob that
decides how badly. Turning it off falls back to the residual count for
the denominator degrees of freedom, which is far cheaper and barely
differs once a cell type has thousands of cells.

`use_cell_quality` conditions on the cell's own quality covariate, which
stage one has already regressed out of the scores by ordinary least
squares. The default conditions on it twice, because upstream does.

## References

Jerby-Arnon & Regev, Nature Biotechnology, 2022
