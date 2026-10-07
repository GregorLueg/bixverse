# CellPhoneDB ligand/receptor analysis

**\[experimental\]** Mirrors the CellPhoneDB v5 numerics. A complex
takes the minimum over its subunits for both the mean and the fraction
expressing. With `statistical = TRUE`, cluster labels are permuted
globally across the included cells and `p = #(perm > real) / n_perm`.
With `deg_genes`, the expression gate is additionally restricted to
interactions where partner a is differentially expressed in the sender
or partner b in the receiver. Reads the normalised layer of the L/R
genes only, once.

## Usage

``` r
rs_sc_cellphonedb(
  f_path_gene,
  partner_a,
  partner_b,
  clusters,
  pair_a,
  pair_b,
  deg_genes,
  statistical,
  params,
  verbose
)
```

## Arguments

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- partner_a:

  List of integer vectors (0-indexed!). Subunit gene indices of partner
  a per interaction; a single gene is a length one vector.

- partner_b:

  List of integer vectors (0-indexed!). Same for partner b.

- clusters:

  List of integer vectors (0-indexed!). Disjoint cell indices per
  cluster. Cells in no cluster are ignored.

- pair_a, pair_b:

  Optional integer vectors (0-indexed!) of equal length giving the
  ordered cluster pairs `(pair_a[i], pair_b[i])` to test. `NULL` tests
  all ordered pairs.

- deg_genes:

  Optional list of integer vectors (0-indexed!), one per cluster in the
  order of `clusters`, with the differentially expressed genes. Replaces
  the expression gate with the DEG gate.

- statistical:

  Boolean. Shall the permutation p-values be computed.

- params:

  List. See
  [`params_sc_cellphonedb()`](https://gregorlueg.github.io/bixverse/reference/params_sc_cellphonedb.md).

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with:

- means - Numeric matrix (interactions x pairs) of interaction means.

- pvals - Numeric matrix (interactions x pairs) of permutation p-values,
  1 where the mean is zero or the gate fails. `NULL` unless
  `statistical = TRUE`.

- gate - Logical matrix (interactions x pairs). The expression gate, or
  the DEG gate if `deg_genes` was supplied.

- pair_a, pair_b - The tested cluster pairs (0-indexed), the columns of
  the matrices above.

- genes - Unique L/R gene indices (0-indexed), the rows of `gene_mean`
  and `gene_pct`.

- gene_mean, gene_pct - Numeric matrices (genes x clusters) with the
  mean normalised expression and the fraction of expressing cells.
