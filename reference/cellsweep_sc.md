# Remove ambient and bulk contamination with CellSweep

Fits the CellSweep multinomial mixture model and writes a new
`SingleCells` object holding the denoised counts. Every observed count
is split three ways by EM: ambient contamination drawn from a
per-emulsion profile, bulk contamination drawn from a global profile,
and true expression drawn from the barcode's cell-type profile.
Subtracting the first two gives the denoised matrix.

## Usage

``` r
cellsweep_sc(
  target,
  input,
  celltype_column,
  sample_column,
  empty_params,
  cellsweep_params = params_sc_cellsweep(),
  sc_qc_param = params_sc_min_quality(),
  csc_mem_gb = NULL,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
  .verbose = TRUE
)
```

## Arguments

- target:

  `SingleCells`. A fresh object pointing at the directory the denoised
  counts should be written to.

- input:

  `SingleCells`. The raw object, with the empty droplets still in it.

- celltype_column:

  String. Obs column with the cell-type labels.

- sample_column:

  String. Obs column identifying the emulsion. One independent fit per
  level.

- empty_params:

  List. See
  [`params_sc_empty_droplets()`](https://gregorlueg.github.io/bixverse/reference/params_sc_empty_droplets.md).
  Required: there is no safe default, since the recommended
  `method = "supplied"` needs the name of the obs column holding the
  mask.

- cellsweep_params:

  List. See
  [`params_sc_cellsweep()`](https://gregorlueg.github.io/bixverse/reference/params_sc_cellsweep.md).

- sc_qc_param:

  List. See
  [`params_sc_min_quality()`](https://gregorlueg.github.io/bixverse/reference/params_sc_min_quality.md).
  Only `target_size` is used, to scale the new normalised layer.

- csc_mem_gb:

  Optional numeric. Memory in GB for the buffers of the cell-to-gene
  (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL` (default)
  converts in one pass and holds the whole matrix. Set a cap for large
  data sets; every extra phase re-reads the cell file once.

- streaming, batch_size, max_genes_in_memory, cell_batch_size:

  Replaced by `csc_mem_gb` and ignored. **\[deprecated\]**

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

The `target` object, populated with the denoised counts. The fit lands
in obs as `cellsweep_alpha`, `cellsweep_z`, `cellsweep_beta`,
`cellsweep_ll` and `cellsweep_converged`, and in var as
`cellsweep_ambient`.

## Details

Two things about where this sits in the workflow, because both are easy
to get wrong.

**It runs after annotation, not before.** The model subtracts against
cell-type profiles, so it needs the labels on input. The chain is:
ingest the raw barcodes, mask and cluster and annotate as usual, then
`cellsweep_sc()`, then redo feature selection and reduction on the clean
counts.

**It needs the empty droplets.** The ambient profile is estimated from
their pooled counts, and they stay in the EM with their contamination
fraction pinned at 1. That means the object has to have been ingested
permissively:

    obj <- load_mtx(
      obj,
      sc_mtx_io_param = get_cell_ranger_params(path),
      sc_qc_param = params_sc_min_quality(
        min_unique_genes = 0L, min_lib_size = 0L, min_cells = 0L
      )
    )

The load-time cutoffs are irreversible, so the defaults
(`min_lib_size = 250L`) delete exactly the barcodes CellSweep trains on.
This function checks for that and errors rather than fitting a garbage
profile.

One EM fit per sample, since the ambient profile is a property of a
single emulsion. `sample_column` is required for that reason; pooling
samples into one ambient profile is wrong even when it runs.

Barcodes partition three ways: empty droplets, annotated barcodes that
passed QC, and everything else. The third group is excluded from the fit
and does not appear in the output.

The denoised counts land in a new directory. The raw layer takes the
stochastically rounded values so the negative binomial methods
downstream still see integers, and the normalised layer keeps the float
magnitudes. Entries whose denoised value rounds to zero are dropped from
both layers: the two share one index set, so keeping them would put
explicit zeros in the raw layer and `library_size` would stop being the
sum of what is stored, breaking every consumer that computes a fraction
of the library.

## References

Sullivan et al., CellSweep, 2025
