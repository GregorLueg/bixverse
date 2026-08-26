# Generate an LdaKSweepResult instance

Takes the raw Rust output of
[`rs_lda_k_sweep()`](https://gregorlueg.github.io/bixverse/reference/rs_lda_k_sweep.md)
and turns the per-`k` metrics into a data.table, keeping the fitted
models as an attribute so the winning one can be pulled out without
refitting.

## Usage

``` r
new_lda_k_sweep_result(sweep_res, doc_ids, term_ids, params)
```

## Arguments

- sweep_res:

  List. The raw return of
  [`rs_lda_k_sweep()`](https://gregorlueg.github.io/bixverse/reference/rs_lda_k_sweep.md).

- doc_ids:

  Character vector. The document identifiers.

- term_ids:

  Character vector. The term identifiers.

- params:

  List. The parameters the sweep was run with.

## Value

An object of class `LdaKSweepResult`, which is also a data.table.
`combined_score` is `NA` for any `k` the coherence topic-count floor
excluded from selection.

## References

Arun, et al., PAKDD, 2010; Cao, et al., Neurocomputing, 2009; Mimno, et
al., EMNLP, 2011
