# Principal component regression on batch

**\[experimental\]** Regresses each embedding dimension on the batch
labels (one-way ANOVA) and weights the per-dimension R-squared by the
variance each dimension carries. Compare the value on the uncorrected
PCA with the one on the corrected embedding, `(pre - post) / pre`,
rather than reading it on its own.

## Usage

``` r
rs_pcr(embedding, batch_vector)
```

## Arguments

- embedding:

  Numeric matrix. Cells x dimensions, ideally a PCA.

- batch_vector:

  Integer vector. The batch per cell. The codes need not be 0-based or
  contiguous.

## Value

A list with the following items

- var_explained - Variance per embedding dimension.

- r_squared - R-squared of batch per embedding dimension.

- pcr - Variance-weighted R-squared of batch.

## References

Büttner, et al., Nat Methods, 2019; Luecken, et al., Nat Methods, 2022
