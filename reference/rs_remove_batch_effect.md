# Remove batch effects from a log-expression matrix

**\[experimental\]** limma's `removeBatchEffect`, via the `edge-rs`
crate. Batch gets sum-to-zero contrasts, is fitted jointly with the
design of interest and only the batch part is subtracted. Meant for
plotting and unsupervised work; for testing put batch into the design.

## Usage

``` r
rs_remove_batch_effect(x, batch, design)
```

## Arguments

- x:

  Numeric matrix. Log-expression values of genes x samples.

- batch:

  Integer vector. Batch label per sample, e.g. `as.integer(factor(x))`.
  Needs at least two distinct labels.

- design:

  Numeric matrix or NULL. The design of interest, samples x
  coefficients, whose effects are protected. NULL is an intercept only.

## Value

Numeric matrix of corrected values, genes x samples.

## References

Smyth, Stat Appl Genet Mol Biol, 2004
