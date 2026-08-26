# Run DIALOGUE over meta cells

**\[experimental\]** The meta cell entry point into DIALOGUE, see
[`rs_dialogue_sc()`](https://gregorlueg.github.io/bixverse/reference/rs_dialogue_sc.md)
for what the method does. This is a shim rather than a second
implementation: everything DIALOGUE asks of the expression matrix is
per-gene, so the in-memory matrix is wrapped as a gene-major reader and
the same core runs.

Only the normalised layer is ever read, so `sparse_data` has to carry
the normalised counts. The `data` layer is cast to integers on the way
in and then goes unused.

Meta cells are already aggregates, so the sample a meta cell belongs to
has to be unambiguous: build them within samples, not across them. The
random intercept in stage two is over samples, and a meta cell
straddling two of them has no well-defined level.

## Usage

``` r
rs_mc_dialogue(
  sparse_data,
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

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`,
  `cs_type`, `nrow` and `ncol`, holding the *normalised* meta cell
  counts with shape (metacells, genes).

- cell_type_indices:

  List of integer vectors. 0-indexed(!) positions of the meta cells
  belonging to each cell type. At least two cell types are needed.

- features:

  List of numeric matrices, one per cell type, shaped
  `n_metacells_in_type x k_i` with rows aligned to `cell_type_indices`.
  Needs at least two columns per cell type.

- sample_ids:

  Integer vector. 0-indexed(!) sample code per meta cell, over all meta
  cells rather than per cell type.

- cell_quality:

  Numeric vector. Quality covariate per meta cell, indexed the same way
  as `sample_ids`.

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

A list, identical in shape to
[`rs_dialogue_sc()`](https://gregorlueg.github.io/bixverse/reference/rs_dialogue_sc.md).

## References

Jerby-Arnon & Regev, Nature Biotechnology, 2022
