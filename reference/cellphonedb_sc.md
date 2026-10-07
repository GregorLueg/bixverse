# CellPhoneDB ligand-receptor analysis

Runs the CellPhoneDB v5 analysis on a `SingleCells` object. For every
interaction and every ordered (sender, receiver) cluster pair, partner a
is measured in the sender and partner b in the receiver. The interaction
mean is the average of the two partner means, zero if either is not
expressed. A complex takes the minimum over its subunits. The gate
requires both partners to be expressed in more than `threshold` of their
cluster's cells.

Three methods, as in CellPhoneDB:

- `"statistical"` - Permutes the cluster labels and reports
  `p = #(perm > real) / n_perm`, 1 where the mean is zero or the gate
  fails.

- `"degs"` - No permutations. The gate also needs partner a to be
  differentially expressed in the sender or partner b in the receiver,
  as given by `deg_table`.

- `"simple"` - Means and the expression gate only.

Reads the normalised layer of the L/R genes only, once. Interactions
with any subunit missing from the object are dropped with a warning.

## Usage

``` r
cellphonedb_sc(
  object,
  celltype_colname,
  lr_db,
  method = c("statistical", "degs", "simple"),
  deg_table = NULL,
  senders = NULL,
  receivers = NULL,
  gene_id_col = "gene_id",
  params = params_sc_cellphonedb(),
  .verbose = TRUE
)
```

## Arguments

- object:

  A `SingleCells` object.

- celltype_colname:

  String. Name of the cluster column in `obs`. Cells with a missing
  label are ignored.

- lr_db:

  data.table. The ligand-receptor database, see
  [`get_cellphonedb_db()`](https://gregorlueg.github.io/bixverse/reference/get_cellphonedb_db.md).
  Needs `interaction_id`, `partner_a`, `partner_b` and the subunit list
  columns `genes_a` and `genes_b`.

- method:

  String. One of `c("statistical", "degs", "simple")`.

- deg_table:

  Optional data.table with the columns `cluster_id` and `gene`, holding
  the differentially expressed genes per cluster, e.g. a filtered
  [`find_all_markers_sc()`](https://gregorlueg.github.io/bixverse/reference/find_all_markers_sc.md)
  result. Required for `method = "degs"`.

- senders, receivers:

  Optional character vectors. Restrict the tested pairs to these sender
  and receiver clusters. `NULL` uses all clusters.

- gene_id_col:

  String. The var column holding the gene symbols used in `lr_db`. A
  duplicated symbol resolves to its first gene. Defaults to `"gene_id"`.

- params:

  List. See
  [`params_sc_cellphonedb()`](https://gregorlueg.github.io/bixverse/reference/params_sc_cellphonedb.md).

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

A data.table with one row per interaction and cluster pair:

- interaction_id, partner_a, partner_b - From `lr_db`.

- sender, receiver - The cluster pair.

- mean - The interaction mean.

- pval - Permutation p-value, `NA` unless `method = "statistical"`.

- gate - The expression gate, or the DEG gate for `method = "degs"`.

## References

Efremova et al., Nat Protoc, 2020; Troulé et al., Nat Protoc, 2025.
