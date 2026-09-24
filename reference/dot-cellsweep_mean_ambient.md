# Average the per-sample ambient profiles

The ambient profile is per-emulsion, so a single var column can only
carry a summary. The unweighted mean across samples is that summary; the
per-sample profiles are not persisted.

## Usage

``` r
.cellsweep_mean_ambient(fits)
```

## Arguments

- fits:

  List. One entry per sample, each with an `ambient` vector.

## Value

Numeric vector, one entry per gene.
