# Wrapper function for the ATAC peak calling parameters

Parameters of
[`call_peaks_sc()`](https://gregorlueg.github.io/bixverse/reference/call_peaks_sc.md),
a native port of MACS3 `callpeak --nomodel` without a control track.

## Usage

``` r
params_atac_call_peaks(
  tag_source = c("insertion", "fragment_start"),
  keep_dup = NULL,
  shift = -100L,
  extsize = 200L,
  genome_size = 2.7e+09,
  llocal = 10000L,
  q_value = 0.05,
  min_length = NULL,
  max_gap = NULL
)
```

## Arguments

- tag_source:

  String. `"insertion"` takes both Tn5 insertions of a fragment,
  `"fragment_start"` the fragment start only, which is what MACS reads
  from the fragment file Signac hands it. One of
  `c("insertion", "fragment_start")`. Defaults to `"insertion"`.

- keep_dup:

  Integer or `NULL`. Tags kept per position. `NULL` keeps all. MACS's
  default, and so Signac's, is 1. Defaults to `NULL`.

- shift:

  Integer. Tag shift before extension (`--shift`). Defaults to `-100L`.

- extsize:

  Integer. Tag extension, the fragment size (`--extsize`). Defaults to
  `200L`.

- genome_size:

  Numeric. Effective genome size (`-g`). The default is Signac's,
  MACS3's own `hs` is 2.913e9. Defaults to `2700000000.0`.

- llocal:

  Integer. Window of the local lambda (`--llocal`). Defaults to
  `10000L`.

- q_value:

  Numeric. q-value cutoff (`-q`). Defaults to `0.05`.

- min_length:

  Integer or `NULL`. Minimum peak length. `NULL` uses `extsize`.
  Defaults to `NULL`.

- max_gap:

  Integer or `NULL`. Largest gap merged into a peak. `NULL` uses MACS's
  tag size: 1 for insertions, the mean length of the group's first
  fragments for fragment starts. Defaults to `NULL`.

## Value

A named list with the following elements:

- tag_source - String. `"insertion"` takes both Tn5 insertions of a
  fragment, `"fragment_start"` the fragment start only, which is what
  MACS reads from the fragment file Signac hands it. One of
  `c("insertion", "fragment_start")`. Defaults to `"insertion"`.

- keep_dup - Integer or `NULL`. Tags kept per position. `NULL` keeps
  all. MACS's default, and so Signac's, is 1. Defaults to `NULL`.

- shift - Integer. Tag shift before extension (`--shift`). Defaults to
  `-100L`.

- extsize - Integer. Tag extension, the fragment size (`--extsize`).
  Defaults to `200L`.

- genome_size - Numeric. Effective genome size (`-g`). The default is
  Signac's, MACS3's own `hs` is 2.913e9. Defaults to `2700000000.0`.

- llocal - Integer. Window of the local lambda (`--llocal`). Defaults to
  `10000L`.

- q_value - Numeric. q-value cutoff (`-q`). Defaults to `0.05`.

- min_length - Integer or `NULL`. Minimum peak length. `NULL` uses
  `extsize`. Defaults to `NULL`.

- max_gap - Integer or `NULL`. Largest gap merged into a peak. `NULL`
  uses MACS's tag size: 1 for insertions, the mean length of the group's
  first fragments for fragment starts. Defaults to `NULL`.

## References

Zhang, et al., Genome Biol, 2008.
