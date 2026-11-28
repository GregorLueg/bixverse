# Wrapper function for the gene activity parameters

Parameters of
[`gene_activity_sc()`](https://gregorlueg.github.io/bixverse/reference/gene_activity_sc.md),
Signac's `GeneActivity()`: per gene, the longest span of its biotype,
under `max_width`, extended `upstream` of the TSS.

## Usage

``` r
params_atac_gene_activity(
  upstream = 2000L,
  downstream = 0L,
  max_width = 500000L,
  biotypes = "protein_coding"
)
```

## Arguments

- upstream:

  Integer. Bases added upstream of the TSS. Defaults to `2000L`.

- downstream:

  Integer. Bases added downstream of the gene end. Defaults to `0L`.

- max_width:

  Integer. Genes this long or longer are dropped. Defaults to `500000L`.

- biotypes:

  String vector or `NULL`. Gene biotypes to keep. `NULL` keeps all.
  Defaults to `"protein_coding"`.

## Value

A named list with the following elements:

- upstream - Integer. Bases added upstream of the TSS. Defaults to
  `2000L`.

- downstream - Integer. Bases added downstream of the gene end. Defaults
  to `0L`.

- max_width - Integer. Genes this long or longer are dropped. Defaults
  to `500000L`.

- biotypes - String vector or `NULL`. Gene biotypes to keep. `NULL`
  keeps all. Defaults to `"protein_coding"`.

## References

Stuart, et al., Nat Methods, 2021.
