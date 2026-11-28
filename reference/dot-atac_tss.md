# TSSs of the GTF genes of the biotypes

TSSs of the GTF genes of the biotypes

## Usage

``` r
.atac_tss(gtf, biotypes)
```

## Arguments

- gtf:

  String. The GTF.

- biotypes:

  Optional character vector. `NULL` keeps all.

## Value

The Rust result: `chrom`, `pos`, `strand`, `gene_id`, `gene_name`.
