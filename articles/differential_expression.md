# Differential expression with donor structure

## Intro

You have 3,000 monocytes from eight donors and you want to know which
genes respond to a treatment. The obvious move is a Wilcoxon test over
cells. It is also the wrong one, and the reason is worth being precise
about.

Those 3,000 cells are not 3,000 independent observations. They came from
eight people. Cells from the same donor are correlated, because donors
differ in genotype, in batch, in how the sample was handled, in
everything. A test that treats each cell as independent has an effective
sample size of eight and thinks it has 3,000, so its p-values are far
too small. This is pseudoreplication, and in single cell data it is the
default failure mode rather than an exotic one.

There are two honest ways out.

**Pseudobulk.** Sum the counts per donor per condition and run a bulk
method. The unit of replication becomes the donor, which is correct. You
throw away everything the cell-level data told you about variability
within a donor, and with few donors you end up with very little power.

**A mixed model.** Keep the cells, but model the donor as a random
effect. [NEBULA](https://doi.org/10.1038/s42003-021-02146-6) fits a
negative binomial gamma mixed model per gene, splitting the variance
into a subject-level term and a cell-level term. You keep the cell-level
information and still pay the right price for the donor structure.

`bixverse` gives you all three, and this vignette runs them side by
side, first on a real effect and then on no effect at all. The second
one is the interesting bit.

``` r

library(bixverse)
library(data.table)
library(SingleCellExperiment)
```

## The data

[Kang, et al.](https://doi.org/10.1038/nbt.4042) pooled PBMCs from eight
lupus donors, split each sample in two, and stimulated one half with
IFN-β. Every donor contributes to both arms, so the design is paired and
the donor is a genuine random effect. Roughly 29k cells before quality
control.

The object is a Bioconductor `SingleCellExperiment`, which
[`load_sce()`](https://gregorlueg.github.io/bixverse/reference/load_sce.html)
takes directly: `colData` becomes obs, `rowData` becomes var, and the
counts go through the usual Rust quality control and normalisation.

Drop the doublets first. The `multiplets` column is the demultiplexing
call, and a doublet is two cells’ counts added together, which is not
something any of these models should see.

``` r

sce <- qs2::qs_read(download_kang_pbmc())
sce <- sce[, sce$multiplets == "singlet" & !is.na(sce$cell)]

dir_kang <- file.path(tempdir(), "kang")
dir.create(dir_kang, showWarnings = FALSE, recursive = TRUE)

sc_object <- SingleCells(dir_data = dir_kang)
sc_object <- load_sce(
  object = sc_object,
  sce = sce,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 200L,
    min_lib_size = 500L,
    min_cells = 20L,
    target_size = 1e4
  ),
  .verbose = TRUE
)
#> Pulling the raw counts out of the SingleCellExperiment.
#> Pulling the obs and var data out of the object
#> Writing counts to disk.
#> Generating gene-based data.
#>  Using light streaming for the CSR to CSC conversion.
#> Writing to the DuckDB.
#> Setting internal mapping.

sc_object
#> Single cell experiment (Single Cells).
#>   No cells (original): 24562
#>    To keep n: 24562
#>   No genes: 12132
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
#>   MAGIC imputed: none
#>   Stale artefacts: none
```

Look at the design before anything else.

``` r

obs <- sc_object[[]]

table(obs$ind, obs$stim)
#>       
#>        ctrl stim
#>   101   858 1166
#>   107   517  525
#>   1015 2738 2352
#>   1016 1723 1635
#>   1039  408  584
#>   1244 1879 1464
#>   1256 2108 2026
#>   1488 2030 2549
```

Eight donors, both arms, hundreds to thousands of cells each. That is
what NEBULA needs. It also needs enough cells per donor: below thirty it
quietly downgrades from its LN method to HL, as the original package
does.

## One cell type at a time

Differential expression across a mixed population mostly measures
composition. If the treatment shifts the monocyte fraction, every
monocyte gene comes out “differential” whether or not any monocyte
changed. Work within a cell type.

``` r

mono <- SingleCellsSubset(
  sc_object,
  grouping_column = "cell",
  group = "CD14+ Monocytes"
)

mono <- find_hvg_sc(mono, hvg_no = 500L, .verbose = FALSE)
hvgs <- get_gene_names_from_idx(mono, get_hvg(mono), rust_based = TRUE)

mono
#> Single cell experiment (subset).
#>   No cells: 5355
#>   No genes: 12132
#>   Group: cell = CD14+ Monocytes
#>   HVG calculated: TRUE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
#>   Stale artefacts: none
```

Restricting to the highly variable genes is not just for speed. NEBULA
fits an optimiser per gene, and running it over 12,000 genes to find the
300 that matter is a poor trade.

## The real contrast

IFN-β stimulation is about as strong a perturbation as single cell data
offers. All three approaches should find it. Same gene set for all
three, so the counts are comparable.

``` r

mono_obs <- mono[[]]

wilcox_res <- find_markers_sc(
  mono,
  cells_1 = mono_obs[stim == "stim", cell_id],
  cells_2 = mono_obs[stim == "ctrl", cell_id],
  .verbose = FALSE
)

wilcox_res[gene_id %in% hvgs & fdr <= 0.05, .N]
#> [1] 146
```

For pseudobulk, the samples are donor by condition, so sixteen of them.
The design blocks on the donor, which is what makes it a paired test.

``` r

cell_list <- split(
  mono_obs$cell_id,
  paste(mono_obs$ind, mono_obs$stim, sep = "_")
)

samples <- data.table(sample = names(cell_list))
samples[, c("ind", "stim") := tstrsplit(sample, "_", fixed = TRUE)]
design_pb <- stats::model.matrix(~ ind + stim, data = samples)

pb_res <- pseudobulk_dge_sc(
  object = mono,
  cell_list = cell_list,
  design = design_pb,
  coef = "stimstim",
  edger_params = params_edger_ql(filter = FALSE),
  .verbose = FALSE
)

pb_res[feature_id %in% hvgs & fdr <= 0.05, .N]
#> [1] 253
```

NEBULA takes the design as a formula against the obs table and the donor
column separately. The donor is not a term in the design; it is the
random effect.

``` r

nebula_res <- nebula_sc(
  object = mono,
  subject_col = "ind",
  design = ~stim,
  genes_to_use = hvgs,
  .verbose = TRUE
)

nebula_res
#> ScNebula: 429 genes fitted, 2 coefficients
#>   Method:   ln | subject column: ind
#>   Tested:   stimstim
#>   Subjects: 8 | cells: 5355
#>   Warnings: 0 did not converge, 89 collapsed to a plain NB
```

``` r

head(
  nebula_res$results[order(p_value), .(gene_id, log_fc, z, p_value, fdr)],
  8
)
#>     gene_id    log_fc         z p_value   fdr
#>      <char>     <num>     <num>   <num> <num>
#> 1: APOBEC3B  4.666642  70.90766       0     0
#> 2:     CCL8  5.968108  84.06402       0     0
#> 3:   CXCL10  5.614680 108.35227       0     0
#> 4:     CTSC  2.641458  46.26615       0     0
#> 5:      IL8 -2.067606 -52.81027       0     0
#> 6:   CXCL11  5.931948  44.58929       0     0
#> 7:   FAM26F  3.230210  54.80266       0     0
#> 8:    IL1RN  4.216905  56.46412       0     0
```

CXCL10, CCL8 and APOBEC3B at the top, IL8 down. That is the interferon
response, and it is reassuring rather than surprising.

Note that the Wilcoxon test finds *fewer* genes here than the other two,
not more. Pseudoreplication inflates significance, but the Wilcoxon test
is also a rank test with a expression-proportion filter, and it is
simply less powerful for a large shift in mean. On a real, strong effect
all three agree the effect is there. The difference shows up somewhere
else.

## No contrast at all

Here is the test that matters. Take the control cells only, so there is
no treatment anywhere in the data, and split the eight donors into two
arbitrary groups of four. Any gene called differential is a false
positive by construction.

``` r

ctrl <- mono_obs[stim == "ctrl"]
donors <- sort(unique(ctrl$ind))
group_a <- donors[c(1, 3, 5, 7)]

ctrl[, fake := ifelse(ind %in% group_a, "a", "b")]

table(ctrl$ind, ctrl$fake)
#>       
#>          a   b
#>   101  202   0
#>   107    0 209
#>   1015 779   0
#>   1016   0 369
#>   1039 112   0
#>   1244   0 415
#>   1256 390   0
#>   1488   0 294
```

``` r

wilcox_null <- find_markers_sc(
  mono,
  cells_1 = ctrl[fake == "a", cell_id],
  cells_2 = ctrl[fake == "b", cell_id],
  .verbose = FALSE
)

wilcox_null[gene_id %in% hvgs & fdr <= 0.05, .N]
#> [1] 45
```

``` r

cell_list_null <- split(ctrl$cell_id, ctrl$ind)

samples_null <- data.table(ind = names(cell_list_null))
samples_null[, fake := ifelse(ind %in% group_a, "a", "b")]
design_null <- stats::model.matrix(~fake, data = samples_null)

pb_null <- pseudobulk_dge_sc(
  object = mono,
  cell_list = cell_list_null,
  design = design_null,
  edger_params = params_edger_ql(filter = FALSE),
  .verbose = FALSE
)

pb_null[feature_id %in% hvgs & fdr <= 0.05, .N]
#> [1] 0
```

[`SingleCellsSubset()`](https://gregorlueg.github.io/bixverse/reference/SingleCellsSubset.md)
subsets a `SingleCells`, not another subset, so fold both restrictions
into one grouping column on the parent.

``` r

parent_obs <- sc_object[[]]
sc_object[["cell_stim"]] <- paste(parent_obs$cell, parent_obs$stim, sep = "__")
sc_object[["fake"]] <- ifelse(parent_obs$ind %in% group_a, "a", "b")

mono_ctrl <- SingleCellsSubset(
  sc_object,
  grouping_column = "cell_stim",
  group = "CD14+ Monocytes__ctrl"
)

nebula_null <- nebula_sc(
  object = mono_ctrl,
  subject_col = "ind",
  design = ~fake,
  genes_to_use = hvgs,
  .verbose = TRUE
)

nebula_null$results[fdr <= 0.05, .N]
#> [1] 6
```

The cell-level test calls dozens of genes differential between two
groups of donors that differ in nothing. The donor-aware tests call
almost none. That gap is the entire argument, and it does not go away
with more cells: sequencing more cells per donor makes the Wilcoxon test
*more* confident about a difference that is not there.

NEBULA is not exactly zero, and it should not be advertised as such.
Eight donors split four against four is a small experiment, real donors
genuinely differ, and an arbitrary split can land partway along some
real axis of donor variation. A handful of calls out of several hundred
genes is the honest expectation.

## Reading the rest of the output

The Wald test is the headline, but the fit carries the variance
decomposition that motivated the model.

``` r

head(
  nebula_res$results[
    order(-subject_overdispersion),
    .(
      gene_id,
      subject_overdispersion,
      cell_overdispersion,
      sigma_at_bound,
      convergence
    )
  ],
  6
)
#>                   gene_id subject_overdispersion cell_overdispersion
#>                    <char>                  <num>               <num>
#> 1:               APOBEC3B              1.7497689           0.5120647
#> 2:               HLA-DRB5              1.3323523           0.6969546
#> 3:                   GJB2              0.9802969           9.5504329
#> 4:               TMEM176B              0.9453621           1.0881471
#> 5: AKR1C1_ENSG00000187134              0.9314292          64.4138787
#> 6:               TMEM176A              0.9152085           2.6871323
#>    sigma_at_bound convergence
#>            <lgcl>       <int>
#> 1:          FALSE           1
#> 2:          FALSE           1
#> 3:          FALSE           1
#> 4:          FALSE           1
#> 5:          FALSE           1
#> 6:          FALSE           1
```

`subject_overdispersion` is how much the gene’s expression moves between
donors, and it is exactly the term a cell-level test pretends is zero.
Genes with a large value are the ones a Wilcoxon test would most
confidently get wrong.

Two columns are diagnostics rather than results. `convergence` at or
below `-20` marks a fit that likely failed, and `sigma_at_bound` marks a
gene where the subject-level variance finished pinned on its lower
bound, meaning the mixed model collapsed to a plain negative binomial.
Neither is fatal, but a sweep where most genes show either is telling
you something about the design.

``` r

nebula_res$results[, .(
  n = .N,
  failed = sum(convergence <= -20L),
  collapsed = sum(sigma_at_bound)
)]
#>        n failed collapsed
#>    <int>  <int>     <int>
#> 1:   429      0        89
```

The fixed effects and their standard errors are on the object as
matrices of genes by coefficients, so the full model is available
without refitting.

``` r

head(nebula_res$coefficients, 4)
#>          (Intercept)  stimstim
#> HBB        -7.310408 0.5106342
#> APOBEC3B   -8.637498 4.6666425
#> HBA2       -8.767250 0.2434404
#> HBA1       -9.472912 0.3307963
```

## Which one to reach for

Use pseudobulk when the number of donors is decent and you want
something conservative, well understood and fast. It is a bulk method
with decades of scrutiny behind it, and
[`run_edger_ql()`](https://gregorlueg.github.io/bixverse/reference/run_edger_ql.md)
gives you the edgeR quasi-likelihood chain in Rust.

Use NEBULA when the cell-level information is worth keeping: few donors,
uneven cell counts, or when you care about the variance decomposition
rather than only the fold change.

Use the Wilcoxon test for marker detection, which is what it is good at.
Asking “which genes distinguish this cluster” is a description of the
data in front of you, not an inference about a population of donors, so
pseudoreplication does not apply. Asking “does this treatment change
expression” is an inference, and then it does.

## Cost

NEBULA fits an optimiser per gene and the fits are independent, so it
parallelises well and scales linearly in genes. On this data 500 genes
over 5,000 cells takes a couple of seconds. Something with 50,000 cells
and 5,000 genes will take minutes rather than seconds, which is fine for
a considered analysis and painful in a loop over twenty cell types.
Restrict the genes.

`gene_batch_size` in
[`params_nebula()`](https://gregorlueg.github.io/bixverse/reference/params_nebula.html)
bounds how much of the store is resident at once. It changes nothing
about the answer, since NEBULA is gene-independent, so leave it alone
unless memory is tight.

## Where next

The differential *abundance* question is the other half of this, and it
has its own
[vignette](https://gregorlueg.github.io/bixverse/articles/differential_abundance.html):
whether the treatment changed which cells are there, rather than what
those cells are expressing. The same donor structure argument applies,
and Milo handles it the same way.
