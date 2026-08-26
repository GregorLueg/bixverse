# Single cell test data with a planted multicellular programme

Generates synthetic data DIALOGUE should be able to solve. Every cell
type gets its own noise and its own sample-level nuisance factors; only
the first feature column and the planted genes carry a shared per-sample
latent, so anything found beyond that is spurious.

## Usage

``` r
generate_dialogue_test_data(
  syn_data_params = params_sc_synthetic_dialogue(),
  seed = 42L
)
```

## Arguments

- syn_data_params:

  List. Contains the parameters for the generation of the synthetic
  data, see:
  [`params_sc_synthetic_dialogue()`](https://gregorlueg.github.io/bixverse/reference/params_sc_synthetic_dialogue.md).

- seed:

  Integer. The seed for the generation of the synthetic data.

## Value

List with the following items

- counts - dgRMatrix with cells x genes of raw counts.

- obs - data.table with `cell_id`, `cell_grp`, `sample_id` and
  `cell_quality`. The last is pure noise, independent of the planted
  programme: hand it to `dialogue_sc(quality_col = "cell_quality")` so
  the covariate does not carry the signal you are trying to find.

- var - data.table with the gene information.

- features - Named list of numeric matrices, one per cell type, with the
  cell identifiers as row names. Hand this straight to
  [`dialogue_sc()`](https://gregorlueg.github.io/bixverse/reference/dialogue_sc.md).

- latent - Numeric vector. The per-sample latent the planted programme
  follows, named by sample.

- planted - Named list of character vectors. The planted gene
  identifiers per cell type.

## Details

Same generator the Rust integration tests use, so the planted structure
cannot drift between the two suites. What differs is the scale: R takes
the raw count layer and lets the normal ingestion path log-normalise it,
where the Rust tests feed the planted layer straight in.

Cells are laid out contiguously by cell type and, within a cell type, by
sample, so every cell type sees every sample. That is the easy case for
the method and is what a fixture should be.

## References

Jerby-Arnon & Regev, Nature Biotechnology, 2022
