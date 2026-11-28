# Link peaks to nearby genes by correlation over groups of cells

ArchR's `addPeak2GeneLinks()`. Peak and gene counts are summed per level
of `group_by` (fine clusters, or meta cell assignments, of the kept
cells), log-normalised per group, `log2(x / total * 1e4 + 1)`, and
correlated over the groups for every peak within `max_distance` of a
gene's TSS. p-values are the two-sided t-test of each correlation, FDR
Benjamini-Hochberg over all pairs. Unlike Signac's `LinkPeaks()` there
is no GC-matched background, which needs the genome sequence.

## Usage

``` r
peak_to_gene_sc(
  object,
  gtf,
  group_by,
  max_distance = 250000L,
  biotypes = "protein_coding",
  genes = NULL,
  min_cells = 10L,
  fdr = 0.05,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCellsMultiModal` class with RNA and ATAC.

- gtf:

  String. Gene annotation for the TSSs, matched to the RNA on `gene_id`
  without its version.

- group_by:

  String. obs column whose levels are the groups; ArchR uses a few
  hundred groups of ~100 cells.

- max_distance:

  Integer. Largest distance from a peak to a TSS.

- biotypes:

  Optional character vector. Gene biotypes to link. `NULL` keeps all.

- genes:

  Optional character vector. Gene ids to restrict to; `NULL` for every
  gene with a TSS.

- min_cells:

  Integer. Groups with fewer kept cells are left out.

- fdr:

  Optional numeric. Keep links below this FDR; `NULL` keeps every tested
  pair.

- .verbose:

  Boolean. Controls verbosity.

## Value

A data.table with `peak_id`, `gene_id`, `gene_name`, `distance` (peak
centre minus TSS), `r`, `p` and `fdr`.

## References

Granja, et al., Nat Genet, 2021
