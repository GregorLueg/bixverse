# Call peaks on the kept cells, pooled or per group

A native port of MACS3 `callpeak --nomodel` without a control track, one
call per group, as Signac's `CallPeaks(group.by = ...)` runs MACS.
Against MACS3 3.1.0 the peak coordinates, summits and pileups are exact
and the scores agree to MACS's printed precision. Pass the result as
`peaks` to
[`add_atac_fragments_sc()`](https://gregorlueg.github.io/bixverse/reference/add_atac_fragments_sc.md)
or
[`load_atac_fragments()`](https://gregorlueg.github.io/bixverse/reference/load_atac_fragments.md)
to count against it.

With a tabix index next to the fragment file (`.tbi`, which Cell Ranger
writes), the file is read one chromosome at a time and the sorted tags
wait in a scratch file (4 bytes per tag) between passes. Without an
index every tag is held at once while reading.

## Usage

``` r
call_peaks_sc(
  object,
  fragments,
  group_by = NULL,
  combine = TRUE,
  peak_params = params_atac_call_peaks(),
  barcode_map = NULL,
  max_mem_gb = NULL,
  tmp_dir = NULL
)
```

## Arguments

- object:

  `SingleCells` class; its kept cells are pooled.

- fragments:

  String. The fragment file.

- group_by:

  Optional string. obs column to call per level of (clusters, samples);
  `NULL` pools all kept cells into one group, `"all"`.

- combine:

  Boolean. Merge overlapping peaks across groups into one set, with the
  groups each was called in. `FALSE` returns every group's peaks.

- peak_params:

  List. Output of
  [`params_atac_call_peaks()`](https://gregorlueg.github.io/bixverse/reference/params_atac_call_peaks.md).

- barcode_map:

  Optional. Cell id to fragment barcode, where they differ.

- max_mem_gb:

  Optional numeric. Tag budget of the reading pass; `NULL` reads one
  chromosome at a time.

- tmp_dir:

  Optional string. Directory of the scratch file; `NULL` for the system
  temporary directory.

## Value

A data.table. With `combine`: `chrom`, `start`, `end` (0-based,
half-open) and `peak_called_in`, the comma-separated groups. Without:
one row per group and peak, with `group`, `summit` and the summit's
`pileup`, `pscore` (`-log10 p`), `fold` and `qscore` (`-log10 q`).

## References

Zhang, et al., Genome Biol, 2008
