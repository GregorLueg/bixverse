# Per-cell QC from the fragment file

Adds `atac_n_fragments`, `atac_nucleosome_signal`, `atac_tss_enrichment`
(Signac's fast definition) and `atac_mito_fraction` to obs for every
cell, plus `atac_frip` on a peak store: its library size over the
fragments, in the store's count unit, so with insertion counts the share
of insertions in peaks (two per fragment) and otherwise Signac's
fragments in peaks over fragments.

The nucleosome signal is mononucleosomal over nucleosome-free fragments,
as in Signac: infinite for a cell with mononucleosomal but no
nucleosome-free fragments (a bad cell, which any upper threshold drops),
and 0 rather than Signac's `NaN` for a cell with neither.

## Usage

``` r
fragment_qc_sc(
  object,
  fragments,
  gtf,
  biotypes = "protein_coding",
  mito_chroms = c("chrM", "MT", "chrMT", "M"),
  barcode_map = NULL
)
```

## Arguments

- object:

  `SingleCellsAtac` or `SingleCellsMultiModal` class.

- fragments:

  String. The fragment file.

- gtf:

  String. Gene annotation for the TSSs.

- biotypes:

  Optional character vector. Gene biotypes whose TSSs count (Signac:
  protein coding). `NULL` keeps all.

- mito_chroms:

  Character vector. Names of the mitochondrial genome.

- barcode_map:

  Optional. Cell id to fragment barcode, where they differ.

## Value

The object with the QC columns added to obs.

## References

Stuart, et al., Nat Methods, 2021
