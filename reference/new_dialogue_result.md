# Constructor for the DialogueResult class

Wraps what
[`rs_dialogue_sc()`](https://gregorlueg.github.io/bixverse/reference/rs_dialogue_sc.md)
and
[`rs_mc_dialogue()`](https://gregorlueg.github.io/bixverse/reference/rs_mc_dialogue.md)
return, mapping every index back into cell, gene, cell type and sample
names.

## Usage

``` r
new_dialogue_result(dialogue_res, prepped, gene_names, source_class, params)
```

## Arguments

- dialogue_res:

  List. The raw Rust output.

- prepped:

  List. Output of
  [`.prep_dialogue_inputs()`](https://gregorlueg.github.io/bixverse/reference/dot-prep_dialogue_inputs.md).

- gene_names:

  Character. All gene identifiers, in index order.

- source_class:

  String. Class the result came off.

- params:

  List. Parameters the run used.

## Value

A `DialogueResult` object, a list with the following items

- programmes - data.table. One row per programme and cell type pair,
  with the empirical p-value and the canonical correlation.

- mcp_cell_types - List per programme of the cell types it spans.

- scores - Named list per cell type. Final programme scores, cells x
  programmes.

- cca_scores - The same for stage one's canonical scores.

- refit_fidelity - Matrix, cell types x programmes. Correlation between
  the two. A low value means the gene-level refit drifted away from the
  programme the decomposition found.

- ws - Named list per cell type of the sparse canonical weights.

- kept_features - Named list per cell type of the feature columns that
  survived the ANOVA filter.

- verdicts - data.table. What the meta-analysis concluded per gene.

- signatures - data.table. The permissive and strict gene lists.

- shared_samples - Character. Samples present in every cell type.

- cell_types - Character. The analysed cell types, in order.

- source_class - String.

- params - List.
