# Gene set enrichment

Gene sets are always a **named list of character vectors**. That's the input
format everywhere in this part of the package. Identifiers must match whatever
your data uses, usually Ensembl gene ids.

Pick by what you have:

| You have | Use |
|---|---|
| a gene list, want enriched sets | `gse_hypergeometric()` |
| many gene lists | `gse_hypergeometric_list()` |
| a gene list, want non-redundant GO | `GeneOntologyElim()` + `gse_go_elim_method()` |
| a ranked statistic vector | `calc_fgsea()` |
| a ranked vector, want non-redundant GO | `fgsea_go_elim()` |
| an expression matrix, want per-sample scores | `calc_gsva()`, `calc_ssgsea()`, `calc_singscore()` |
| several contrasts at once | `calc_mitch()` |

## Over-representation

```r
res <- gse_hypergeometric(
  target_genes  = target_genes,
  gene_set_list = h_gene_sets_ls
)
```

The universe defaults to the union of all genes in `gene_set_list`. That's
usually fine, but check it if your assay covers a non-standard gene space,
because a wrong universe quietly shifts every p-value.

Many target sets at once, parallelised in Rust:

```r
res <- gse_hypergeometric_list(
  target_gene_list = list(set_1 = ..., set_2 = ...),
  gene_set_list    = h_gene_sets_ls
)
```

## GO-aware enrichment: the elimination method

Running a hypergeometric per GO term independently gives redundant results,
because significant child terms mechanically inflate their ancestors. The
elimination method walks the DAG from the leaves up and strips the genes of a
significant term out of its ancestors before testing them, so ancestors get
tested on residual signal.

```r
go_data_dt <- get_go_data_human()
go_obj     <- GeneOntologyElim(go_data_dt, min_genes = 3L)

res  <- gse_go_elim_method(go_obj, target_genes = target_genes)
resl <- gse_go_elim_method_list(go_obj, target_gene_list = target_list)
```

Human GO ships with the package as parquet in `inst/extdata`, so this needs no
network. `load_go_human_data()` gives you the raw three tables,
`get_go_data_human()` the processed version, `process_go_data()` if you're
assembling your own.

The DAG traversal is sequential per gene set but parallel within a level. 100
target sets against full human GO takes seconds.

### Or simplify after the fact

The alternative is an unconstrained test over all GO terms, then collapsing
redundancy with semantic similarity:

```r
res_simplified <- simplify_hypergeom_res(
  res             = go_results_unfiltered,
  parent_child_dt = go_parent_child_dt,
  weights         = setNames(c(0.8, 0.6), c("is_a", "part_of"))
)
```

This works on any enrichment result, not just bixverse's, which is the reason to
prefer it. It does not adjust the test statistics the way elimination does,
which is the reason not to.

## Rank-based (GSEA)

```r
res <- calc_fgsea(
  stats       = ranked_vector,     # named numeric, names are genes
  pathways    = pathway_list,
  gsea_params = params_gsea(min_size = 15L)
)
```

`calc_fgsea()` is the multilevel adaptive method, `calc_fgsea_simple()` the
simple permutation version, `calc_gsea_traditional()` the original Subramanian
formulation. All three agree closely with the reference implementations, they
just run a lot faster.

GO-aware versions: `fgsea_go_elim(go_obj, stats = ranked_vector)` and
`fgsea_simple_go_elim()`.

## Per-sample scoring

For an expression matrix (genes in rows, samples in columns) where you want one
score per sample per pathway.

**GSVA**, two kernels. Gaussian for continuous data, Poisson for raw counts:

```r
calc_gsva(exp = X,        pathways = gs, kernel = "gaussian")
calc_gsva(exp = X_counts, pathways = gs, kernel = "poisson")
```

**ssGSEA**:

```r
calc_ssgsea(exp = X, pathways = gs)
```

**singscore**, which ranks once then scores many times:

```r
ranks <- calc_singscore_rank(exp = X)

# one signature, optionally with a permutation null
calc_singscore(ranks, up_set = up_set, down_set = down_set,
               n_permutations = 1000L, seed = 42L)

# many signatures, returns $scores and $dispersions
calc_singscore_multi(ranks, up_pathways = gs[1:10])
```

For a targeted panel (NanoString, RT-qPCR) rather than a transcriptome-wide
assay, hand `calc_singscore_rank()` a set of stable genes:

```r
calc_singscore_rank(exp = X, stable_genes = c("g1", "g2", "g3"))
```

`params_gsva()` and `params_ssgsea()` carry the tuning for the first two.

## Multi-contrast (mitch)

When you have several contrasts and want pathways moving jointly:

```r
res <- calc_mitch(
  contrast_mat  = contrast_data,   # genes x contrasts, numeric
  gene_set_list = gene_sets
)
```

Returns per-pathway MANOVA p-values, FDR, and the contrast-level scores.
Pathways below the internal minimum size are dropped from the output rather than
returned with `NA`.

## Getting gene sets

`msigdbr` is the usual source, and it's in `Suggests`:

```r
h <- msigdbr::msigdbr(species = "human", collection = "H")
gs <- split(h$ensembl_gene, h$gs_name)
```

For anything CisTarget-related (motif rankings for SCENIC) see
`single-cell-analysis.md`.
