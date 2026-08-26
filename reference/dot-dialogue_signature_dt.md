# Flatten the DIALOGUE signature lists into a table

Flatten the DIALOGUE signature lists into a table

## Usage

``` r
.dialogue_signature_dt(signatures, cell_types, gene_names, list_name)
```

## Arguments

- signatures:

  Nested list `[cell_type][programme]` of `up` / `down` 0-indexed gene
  positions, as Rust returns them.

- cell_types:

  Character. Cell type names, in index order.

- gene_names:

  Character. All gene identifiers, in index order.

- list_name:

  String. Which list this is, `"permissive"` or `"strict"`.

## Value

A data.table with `cell_type`, `programme`, `gene_id`, `direction` and
`list`.
