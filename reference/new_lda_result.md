# Generate an LdaResult instance

Takes the raw Rust output of
[`rs_lda()`](https://gregorlueg.github.io/bixverse/reference/rs_lda.md)
and names both matrices. Rust returns the document-topic matrix as
`k x n_documents` because a document's topic vector is one contiguous
column there; this flips it so both matrices are entity-major, matching
the NMF result classes.

## Usage

``` r
new_lda_result(lda_res, doc_ids, term_ids, params)
```

## Arguments

- lda_res:

  List. The raw return of
  [`rs_lda()`](https://gregorlueg.github.io/bixverse/reference/rs_lda.md).

- doc_ids:

  Character vector. The document identifiers, in row order of the input
  matrix.

- term_ids:

  Character vector. The term identifiers, in column order of the input
  matrix.

- params:

  List. The parameters the model was fitted with.

## Value

An object of class `LdaResult`.
