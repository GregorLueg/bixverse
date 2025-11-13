# Helper to process CisTarget results

Helper to process CisTarget results

## Usage

``` r
process_cistarget_res(cs_ls, gs_name, represented_motifs, represented_genes)
```

## Arguments

- cs_ls:

  List. The result list from the Rust wrapper.

- gs_name:

  String. Name of the tested gene set.

- represented_motifs:

  Character vector. The represented motifs in the rankings.

- represented_genes:

  Character vector. The represented genes in the rankings.

## Value

A data.table with the results if there were any significant motifs.
