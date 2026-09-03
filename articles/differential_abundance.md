# Differential abundance with Milo and MELD

## Intro

Differential expression asks whether cells changed what they express.
Differential abundance asks whether the cells that are there changed at
all. Both matter, and the second is the one people usually answer badly.

The obvious approach is to cluster, count cells per cluster per sample,
and test the proportions. It works, and it inherits every problem the
clustering had. If a condition shifts a subpopulation that your
resolution parameter merged into a bigger cluster, you see nothing. If
it shifts half a cluster, the proportion barely moves. The answer
depends on a choice you made before you asked the question.

Two methods here avoid that, in different ways.

[Milo](https://doi.org/10.1038/s41587-021-01033-z) tests overlapping
neighbourhoods of the kNN graph instead of clusters. Each neighbourhood
is a cell and its neighbours, so the resolution is the graph’s rather
than the clustering’s, and a shift that spans part of a cluster shows
up. Counts per neighbourhood per sample go through the same negative
binomial machinery bulk RNA-seq uses.

[MELD](https://doi.org/10.1038/s41587-020-00803-5) does not test
anything. It smooths the condition labels over the graph and returns,
per cell, how likely that cell is under each condition. It is a
continuous score rather than a set of calls, which suits a graded
response better than a discrete one.

They answer the same question and they should agree. At the end of this
vignette they do, which is a better check on both than either provides
alone.

``` r

library(bixverse)
library(data.table)
library(SingleCellExperiment)
```

## The data

[Baran-Gale, et al.](https://doi.org/10.1242/dev.183996) profiled mouse
thymic epithelial cells at one, four and sixteen weeks. The thymus
involutes with age, so the cell type proportions genuinely move, which
is what makes this a sensible differential abundance example.

> **A stimulation experiment is a bad Milo example**
>
> It is tempting to reuse a treatment dataset here. Do not. In vitro
> stimulation moves cells a long way in embedding space without moving
> the cell type proportions much, so stimulated and control cells
> separate into different neighbourhoods and almost every neighbourhood
> comes out significant. The result looks spectacular and says nothing
> beyond “the treatment did something”, which you knew.
>
> Differential abundance needs a design where the *composition*
> plausibly changed.

The object carries no gene names at all, with the identifiers sitting in
the `rowData` instead.
[`load_sce()`](https://gregorlueg.github.io/bixverse/reference/load_sce.html)
falls back to the first metadata column and tells you it did. 336 genes
have no annotation whatsoever and get a generated identifier, which it
also tells you. Both warnings are expected here.

``` r

sce <- qs2::qs_read(download_thymus_ageing())

dir_thymus <- file.path(tempdir(), "thymus")
dir.create(dir_thymus, showWarnings = FALSE, recursive = TRUE)

sc_object <- SingleCells(dir_data = dir_thymus)
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
#> Warning in flatten(SummarizedExperiment::rowData(sce), rownames(sce), id_label
#> = "gene_id", : No gene names on the object. Using 'ensembl_gene_id' from the
#> metadata instead.
#> Warning in flatten(SummarizedExperiment::rowData(sce), rownames(sce), id_label
#> = "gene_id", : 336 gene(s) have no identifier. Generating one for each.
#> Writing counts to disk.
#> Generating gene-based data.
#>  Using light streaming for the CSR to CSC conversion.
#> Writing to the DuckDB.
#> Setting internal mapping.

sc_object
#> Single cell experiment (Single Cells).
#>   No cells (original): 69180
#>    To keep n: 69180
#>   No genes: 22719
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
#>   MAGIC imputed: none
#>   Stale artefacts: none
```

Now the design, and there is a trap in it. `samp_id` looks like the
sample. It is not: it is the sequencing run, and every run carries all
three ages, because the ages were hashed together and demultiplexed.

``` r

obs <- sc_object[[]]

table(obs$samp_id, obs$age)
#>          
#>            Wk1 Wk16  Wk4
#>   1stRun1 5160 3798 3902
#>   1stRun2 3360 4199 2355
#>   2ndRun1 4017 3728 4629
#>   2ndRun2 3897 3630 4551
#>   3rdRun1 3263 2981 3771
#>   3rdRun2 3745 3284 4910
```

The sample is the run and the age together, which gives eighteen. Using
`samp_id` as the sample would ask Milo whether abundance differs between
sequencing runs, and every run contains every age, so the answer would
be a flat no.

``` r

sc_object[["sample"]] <- paste(obs$samp_id, obs$age, sep = "_")
obs <- sc_object[[]]

table(obs$sample, obs$age)[1:6, ]
#>               
#>                 Wk1 Wk16  Wk4
#>   1stRun1_Wk1  5160    0    0
#>   1stRun1_Wk16    0 3798    0
#>   1stRun1_Wk4     0    0 3902
#>   1stRun2_Wk1  3360    0    0
#>   1stRun2_Wk16    0 4199    0
#>   1stRun2_Wk4     0    0 2355
```

## The embedding

Milo and MELD both run on the kNN graph, so the usual pipeline comes
first.

``` r

sc_object <- find_hvg_sc(sc_object, hvg_no = 2000L, .verbose = FALSE)
sc_object <- calculate_pca_sc(sc_object, no_pcs = 30L, .verbose = FALSE)
```

The default `k` of 15 is fine for clustering and too small for Milo. A
neighbourhood holds `k + 1` cells, and those get spread over the
samples. With eighteen samples and `k = 15` that is under one cell per
sample per neighbourhood, and a negative binomial fitted to counts that
are mostly zero or one has very little to work with.

Rough rule: pick `k` so the average neighbourhood holds at least a few
cells per sample. Eighteen samples and `k = 60` gives around three,
which is enough here.

``` r

sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(k = 60L)
  ),
  .verbose = TRUE
)
#> 
#> Generating sNN graph (full: TRUE).
#> Transforming sNN data to igraph.
```

## Milo

Sampling the neighbourhoods and counting the cells per sample is one
call. The counts come back as neighbourhoods by samples, with the
graph-overlap weighting and the k-th neighbour distances taken at the
same time, since all three are functions of the neighbourhood matrix
alone.

``` r

milo_obj <- get_miloR_abundances_sc(
  object = sc_object,
  sample_id_col = "sample",
  miloR_params = params_sc_miloR(prop = 0.1),
  .verbose = FALSE
)

dim(milo_obj$sample_counts)
#> [1] 5183   18
mean(rowSums(milo_obj$sample_counts))
#> [1] 61
```

The test needs a design table with one row per sample, and its rownames
have to cover the sample names on the counts. Blocking on the run is
worth it: every run contributes all three ages, so run-to-run variation
is a nuisance term the model can remove.

``` r

samp_design <- unique(obs[, .(sample, samp_id, age)])

design_df <- data.frame(
  age = samp_design$age,
  run = samp_design$samp_id,
  row.names = samp_design$sample
)
design_df <- design_df[colnames(milo_obj$sample_counts), , drop = FALSE]

table(design_df$age, design_df$run)
#>       
#>        1stRun1 1stRun2 2ndRun1 2ndRun2 3rdRun1 3rdRun2
#>   Wk1        1       1       1       1       1       1
#>   Wk16       1       1       1       1       1       1
#>   Wk4        1       1       1       1       1       1
```

> **Check which coefficient you are testing**
>
> [`test_nhoods()`](https://gregorlueg.github.io/bixverse/reference/test_nhoods.md)
> defaults to the last column of the design, the way edgeR and limma do.
> With `age` as a factor, R orders the levels alphabetically: `Wk1`,
> `Wk16`, `Wk4`. The last coefficient is therefore **Wk4 against Wk1**,
> not the sixteen week contrast you probably wanted.
>
> Name the coefficient rather than trusting the default.

``` r

milo_obj <- test_nhoods(
  milo_obj,
  design = ~ run + age,
  design_df = design_df,
  coef = "ageWk16"
)

da_res <- get_differential_abundance_res(milo_obj)

da_res[SpatialFDR <= 0.1, .N]
#> [1] 2313
```

``` r

head(da_res[order(PValue), .(Nhood, logFC, F, PValue, FDR, SpatialFDR)], 6)
#>    Nhood     logFC        F       PValue          FDR   SpatialFDR
#>    <int>     <num>    <num>        <num>        <num>        <num>
#> 1:  3080  6.227930 48.12281 4.309899e-12 2.233821e-08 1.433629e-08
#> 2:  2474  5.892283 40.13958 2.490099e-10 5.291019e-07 3.611810e-07
#> 3:  1238  4.911830 39.73336 3.062523e-10 5.291019e-07 3.611810e-07
#> 4:  3347  4.885818 37.70221 8.627143e-10 1.117862e-06 8.274642e-07
#> 5:  1406  4.354272 36.44494 1.639166e-09 1.307584e-06 1.122626e-06
#> 6:   228 -5.949431 36.28411 1.779562e-09 1.307584e-06 1.122626e-06
```

A neighbourhood is not a cell type, so on its own a list of
neighbourhood indices is not much use.
[`add_nhoods_info()`](https://gregorlueg.github.io/bixverse/reference/add_nhoods_info.md)
tags each with the cell type its cells mostly belong to.

``` r

milo_obj <- add_nhoods_info(milo_obj, cell_info = obs$cluster_annot)
da_res <- get_differential_abundance_res(milo_obj)

sig <- da_res[SpatialFDR <= 0.1]
sig[, direction := ifelse(logFC > 0, "up at Wk16", "down at Wk16")]

dcast(
  sig[, .N, by = .(majority_celltype, direction)],
  majority_celltype ~ direction,
  value.var = "N",
  fill = 0
)[order(-`up at Wk16`)]
#>      majority_celltype down at Wk16 up at Wk16
#>                 <char>        <int>      <int>
#>  1: Intertypical.TEC.2           25        266
#>  2: Intertypical.TEC.4          106        256
#>  3: Intertypical.TEC.1          489        227
#>  4: Intertypical.TEC.3           28        196
#>  5:             mTEC.2           44         35
#>  6:             cTEC.2           46         34
#>  7:             mTEC.1          198         18
#>  8:             mTEC.6            0         18
#>  9:             mTEC.3            0         15
#> 10:    PostAire.mTEC.1           11          8
#> 11:        Tuft.mTEC.2            9          8
#> 12:             mTEC.4            0          8
#> 13:              eTEC1            0          5
#> 14:          New.TEC.2            0          2
#> 15:        Tuft.mTEC.1            0          1
#> 16:    PostAire.mTEC.2            1          0
#> 17:       Prolif.TEC.2           88          0
#> 18:       Prolif.TEC.3           94          0
#> 19:         Sca1.TEC.1            9          0
#> 20:             cTEC.1           33          0
#> 21:             mTEC.5           32          0
#> 22:             mTEC.7            3          0
#>      majority_celltype down at Wk16 up at Wk16
#>                 <char>        <int>      <int>
```

The two proliferating populations are the most strongly depleted of
anything here, which is the clearest ageing signal the tissue has: an
involuting thymus stops making new cells. The intertypical TECs dominate
the enriched side.

The medullary subsets split both ways, so “mTECs go up” would be the
wrong summary to take away. That is the point of testing neighbourhoods
rather than clusters: a cluster level test would have averaged those
subsets into one number and reported whichever direction happened to
win.

### On the spatial FDR

Neighbourhoods overlap, so their tests are not independent and a plain
Benjamini-Hochberg adjustment is anti-conservative. Milo’s spatial FDR
weights each p-value by the reciprocal of its connectivity and runs the
step-up on those weights.

``` r

data.table(
  plain_fdr = da_res[FDR <= 0.1, .N],
  spatial_fdr = da_res[SpatialFDR <= 0.1, .N]
)
#>    plain_fdr spatial_fdr
#>        <int>       <int>
#> 1:      2290        2313
```

On this data the two barely differ. That is worth knowing rather than
disappointing: the correction bites when neighbourhoods overlap heavily
and the connectivity varies a lot between them, and here it does not. Do
not assume it is always doing work.

The weighting scheme is a choice. `"k-distance"` weights by the distance
to the k-th neighbour, `"graph-overlap"` by how many cells a
neighbourhood shares with the others. Both are precomputed, so switching
is cheap.

## MELD

MELD needs no design and no neighbourhood sampling. Hand it the
condition column and it smooths the indicator over the graph.

``` r

meld_res <- meld_sc(
  object = sc_object,
  sample_id_col = "age",
  .verbose = FALSE
)

dim(meld_res$norm_scores)
#> [1] 69180     3
head(round(meld_res$norm_scores, 3), 4)
#>                                             Wk1  Wk16   Wk4
#> Ageing_ZsG_1stRun1_HTO_AAACCCAAGATAGCAT-1 0.299 0.384 0.317
#> Ageing_ZsG_1stRun1_HTO_AAACCCAAGCTGACCC-1 0.498 0.183 0.319
#> Ageing_ZsG_1stRun1_HTO_AAACCCAAGGGAGAAT-1 0.350 0.348 0.302
#> Ageing_ZsG_1stRun1_HTO_AAACCCACAAGACCGA-1 0.197 0.492 0.311
```

The normalised scores are clamped at zero and L1 normalised per cell, so
each row is a likelihood over the three ages. A cell sitting at 0.6 for
Wk16 lives in a part of the manifold that is enriched for sixteen week
old thymus.

Averaging per cell type turns that into something readable.

``` r

meld_dt <- data.table(
  cluster = obs$cluster_annot,
  as.data.table(meld_res$norm_scores)
)

per_celltype <- meld_dt[, lapply(.SD, mean), by = cluster]

head(per_celltype[order(-Wk16)], 8)
#>               cluster       Wk1      Wk16       Wk4
#>                <char>     <num>     <num>     <num>
#> 1:             mTEC.6 0.1753658 0.5715476 0.2530866
#> 2:             mTEC.4 0.2843784 0.4353154 0.2803062
#> 3: Intertypical.TEC.2 0.2679963 0.4078628 0.3241410
#> 4: Intertypical.TEC.4 0.3059496 0.3969262 0.2971241
#> 5:          New.TEC.2 0.3001987 0.3930484 0.3067529
#> 6:             mTEC.3 0.3243268 0.3699300 0.3057432
#> 7: Intertypical.TEC.3 0.2868250 0.3655374 0.3476376
#> 8:          New.TEC.1 0.3706154 0.3639596 0.2654250
```

## Do they agree

They are different methods on the same graph, so this is a real check
rather than a formality. Milo’s mean log fold change per cell type
against MELD’s mean Wk16 likelihood per cell type.

``` r

milo_per_ct <- sig[,
  .(milo_logfc = mean(logFC)),
  by = .(cluster = majority_celltype)
]
comparison <- merge(per_celltype, milo_per_ct, by = "cluster")

comparison[order(-milo_logfc), .(cluster, Wk1, Wk16, milo_logfc)]
#>                cluster       Wk1      Wk16  milo_logfc
#>                 <char>     <num>     <num>       <num>
#>  1: Intertypical.TEC.2 0.2679963 0.4078628  2.01249328
#>  2:             mTEC.6 0.1753658 0.5715476  2.00854230
#>  3: Intertypical.TEC.3 0.2868250 0.3655374  1.71859831
#>  4:          New.TEC.2 0.3001987 0.3930484  1.58159637
#>  5:             mTEC.3 0.3243268 0.3699300  1.54551383
#>  6:              eTEC1 0.3348530 0.3477823  1.34009830
#>  7:             mTEC.4 0.2843784 0.4353154  1.29976300
#>  8:        Tuft.mTEC.1 0.3341064 0.2995690  1.28990557
#>  9: Intertypical.TEC.4 0.3059496 0.3969262  0.96017081
#> 10:             mTEC.2 0.3794905 0.2838506 -0.08258995
#> 11:        Tuft.mTEC.2 0.3228378 0.2523843 -0.19504211
#> 12:    PostAire.mTEC.1 0.3763538 0.2730967 -0.22075283
#> 13:             cTEC.2 0.3051083 0.3374357 -0.31610500
#> 14: Intertypical.TEC.1 0.3356320 0.3118319 -1.01757530
#> 15:             mTEC.7 0.3865443 0.2801945 -1.37043196
#> 16:    PostAire.mTEC.2 0.3368777 0.3434280 -1.42151750
#> 17:             mTEC.1 0.4472890 0.2134771 -1.62957261
#> 18:         Sca1.TEC.1 0.4230381 0.1958383 -1.63619727
#> 19:             mTEC.5 0.4041701 0.2323975 -1.67055463
#> 20:             cTEC.1 0.4692116 0.2623401 -1.71274280
#> 21:       Prolif.TEC.2 0.4808143 0.1788700 -1.94402964
#> 22:       Prolif.TEC.3 0.5166495 0.1540147 -2.44038283
#>                cluster       Wk1      Wk16  milo_logfc
#>                 <char>     <num>     <num>       <num>
```

``` r

cor(comparison$Wk16, comparison$milo_logfc, method = "spearman")
#> [1] 0.8644833
```

Strong agreement, and the disagreements are informative rather than
embarrassing. Milo tests neighbourhoods and reports only the ones that
clear a threshold, so a cell type that shifted a little everywhere can
score low on the Milo axis while MELD still registers the shift. The two
are measuring related but not identical things.

## Which one to reach for

Milo when you want calls with error control: a list of regions that
changed, with p-values you can defend and a multiple testing correction
that accounts for the overlap. It needs a real design, enough samples to
fit a negative binomial, and a `k` large enough that the counts are not
mostly zero.

MELD when the response is graded, when you want a per-cell score to
carry into downstream analysis, or when the design is too thin for a
test. It gives you no p-values, which is a feature when you have four
samples and would only have been fooling yourself.

Running both costs almost nothing once the graph is built, and their
agreement is the most useful diagnostic either one has.

## Where next

The other half of this is differential *expression* with the same donor
structure, which has its own
[vignette](https://gregorlueg.github.io/bixverse/articles/differential_expression.html).
Neighbourhood level expression testing, rather than abundance, is not
wired up yet.
