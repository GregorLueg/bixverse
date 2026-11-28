# Call peaks per group of cells

**\[experimental\]** A port of MACS3 `callpeak` without control.

## Usage

``` r
rs_atac_call_peaks(
  fragments,
  barcodes,
  groups,
  tag_source,
  keep_dup,
  shift,
  extsize,
  genome_size,
  llocal,
  q_value,
  min_length,
  max_gap,
  max_mem_gb,
  tmp_dir
)
```

## Arguments

- fragments:

  String. The fragment file.

- barcodes:

  Character vector. Cells.

- groups:

  Integer vector. Group of each cell, `0..n_groups`.

- tag_source:

  String. `"insertion"` or `"fragment_start"`.

- keep_dup:

  Optional integer. Tags kept per position, `NULL` for all.

- shift:

  Integer. Tag shift before extension.

- extsize:

  Integer. Tag extension.

- genome_size:

  Numeric. Effective genome size.

- llocal:

  Integer. Local lambda window.

- q_value:

  Numeric. q-value cutoff.

- min_length:

  Optional integer. Minimum peak length, `NULL` for `extsize`.

- max_gap:

  Optional integer. Largest merged gap, `NULL` for MACS's tag size.

- max_mem_gb:

  Optional numeric. Tag budget of the first pass, `NULL` for one
  chromosome at a time.

- tmp_dir:

  Optional string. Directory of the tag scratch file, `NULL` for the
  system one.

## Value

A list with `group` (0-based), `chrom`, `start`, `end`, `summit`
(0-based), `pileup`, `pscore`, `fold` and `qscore`. Groups in order,
peaks by chromosome name and position.
