# Drop the peaks that overlap a blacklist region

Drop the peaks that overlap a blacklist region

## Usage

``` r
remove_blacklisted(peaks, blacklist)
```

## Arguments

- peaks:

  data.frame. Regions with `chrom`, `start`, `end`, e.g. from
  [`call_peaks_sc()`](https://gregorlueg.github.io/bixverse/reference/call_peaks_sc.md).

- blacklist:

  A BED file (plain or gzip) or a data.frame of regions, e.g. the ENCODE
  blacklist.

## Value

The peaks that touch no blacklist region, in order.
