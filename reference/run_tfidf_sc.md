# Switch the ATAC store to TF-IDF over the kept cells

Signac's method 1, `log1p(x / lib_size * scale * idf)`, with the IDF
over the kept cells. Header metadata only; no count is rewritten. Run
again after the kept cells change. Everything in the ATAC cache (LSI,
embeddings, kNN) was built on the old norm and is dropped; rerun
[`calculate_lsi_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_lsi_sc.md).

## Usage

``` r
run_tfidf_sc(object, scale = 10000)
```

## Arguments

- object:

  `SingleCellsAtac` or `SingleCellsMultiModal` class with ATAC.

- scale:

  Numeric. Scale factor of the term frequencies.

## Value

The object with the store switched to TF-IDF.

## References

Stuart, et al., Nat Methods, 2021
