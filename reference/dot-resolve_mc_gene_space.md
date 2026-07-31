# Resolve the target gene space for a meta cell merge

Resolve the target gene space for a meta cell merge

## Usage

``` r
.resolve_mc_gene_space(gene_lists, feature_space)
```

## Arguments

- gene_lists:

  List of character vectors with the gene identifiers.

- feature_space:

  String. One of `c("intersect", "union")`.

## Value

Character vector with the target gene identifiers, ordered by the first
input.
