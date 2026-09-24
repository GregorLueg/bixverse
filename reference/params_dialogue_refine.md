# Wrapper function for the DIALOGUE refinement parameters

Stage three of DIALOGUE: the cross-partner meta-analysis that decides
which genes survive, and the non-negative refit of the programme scores
onto them.

## Usage

``` r
params_dialogue_refine(
  support_p = 0.1,
  min_support_fraction = 0.333333333333333,
  min_stratum = 5L,
  early_stop_cor = 0.95,
  permissive_p = 0.001,
  strict_p = 0.05
)
```

## Arguments

- support_p:

  Numeric. Adjusted p below which one partner counts as supporting a
  gene. Must be in `(0, 1]`. Defaults to `0.1`.

- min_support_fraction:

  Numeric. Minimum supporting fraction for a stratum to enter the staged
  fit. Must be in `[0, 1]`. Defaults to `0.3333333333333333`.

- min_stratum:

  Integer. Minimum genes in a stratum before it is worth fitting.
  Defaults to `5L`.

- early_stop_cor:

  Numeric. Correlation between the original score and the running fit at
  which the staged fit stops early. Must be in `(0, 1]`. Defaults to
  `0.95`.

- permissive_p:

  Numeric. Fisher-combined p for the permissive gene list, where a gene
  is carried by partner support rather than by a positive coefficient.
  Must be in `(0, 1]`. Defaults to `0.001`.

- strict_p:

  Numeric. Fisher-combined p for the strict gene list, which also
  demands that every partner supports the gene. Must be in `(0, 1]`.
  Defaults to `0.05`.

## Value

A named list with the following elements:

- support_p - Numeric. Adjusted p below which one partner counts as
  supporting a gene. Must be in `(0, 1]`. Defaults to `0.1`.

- min_support_fraction - Numeric. Minimum supporting fraction for a
  stratum to enter the staged fit. Must be in `[0, 1]`. Defaults to
  `0.3333333333333333`.

- min_stratum - Integer. Minimum genes in a stratum before it is worth
  fitting. Defaults to `5L`.

- early_stop_cor - Numeric. Correlation between the original score and
  the running fit at which the staged fit stops early. Must be in
  `(0, 1]`. Defaults to `0.95`.

- permissive_p - Numeric. Fisher-combined p for the permissive gene
  list, where a gene is carried by partner support rather than by a
  positive coefficient. Must be in `(0, 1]`. Defaults to `0.001`.

- strict_p - Numeric. Fisher-combined p for the strict gene list, which
  also demands that every partner supports the gene. Must be in
  `(0, 1]`. Defaults to `0.05`.

## Details

Two gene lists come out. The permissive one asks only for a
Fisher-combined p below `permissive_p`; the strict one is looser on the
p-value but also demands that *every* partner supports the gene. They
are not nested by threshold, they are nested by evidence, and the strict
list is the one to quote.

## References

Jerby-Arnon & Regev, Nature Biotechnology, 2022
