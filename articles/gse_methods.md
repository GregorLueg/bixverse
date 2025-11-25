# Gene Set Enrichment Methods

## Gene Set Enrichment Methods in bixverse

This vignette shows you how to use the different implement gene set
enrichment methods in `bixverse`.

``` r
if (!requireNamespace("fgsea", quietly = TRUE)) BiocManager::install("fgsea")
if (!requireNamespace("msigdbr", quietly = TRUE)) install.packages("msigdbr")
library(bixverse)
library(data.table)
library(magrittr)
```

### Hypergeometric tests

`bixverse` has implementations of the standard hypergeometric tests to
identify gene set enrichment for a given set of target genes. In this
example, we will be loading the Hallmark gene sets from `msigdbr`. Let’s
first load these in and create a list of gene sets (the expected format
of most functions in `bixverse` for these type of analysis.)

``` r
h_gene_sets <- msigdbr::msigdbr(species = "human", collection = "H")

h_gene_sets_ls <- split(h_gene_sets$ensembl_gene, h_gene_sets$gs_name)

set.seed(123L)

target_genes_1 <- c(
  sample(h_gene_sets_ls[["HALLMARK_MYC_TARGETS_V1"]], 25),
  sample(h_gene_sets_ls[["HALLMARK_MYC_TARGETS_V2"]], 25)
)
```

To run the hypergeometric test against one target set, you can just use
`gse_hypergeometric`. The function has parameters to define a minimum
FDR threshold and minimum overlap science. Additionally, the gene
universe will be set (if not provided) to all of the genes represented
in the `gene_set_list`. You can provide your own gene universe if you
wish to.

``` r
results <- gse_hypergeometric(
  target_genes = target_genes_1,
  gene_set_list = h_gene_sets_ls
)

head(results)
#>              gene_set_name odds_ratios        pvals          fdr  hits
#>                     <char>       <num>        <num>        <num> <num>
#> 1: HALLMARK_MYC_TARGETS_V2  170.749267 5.847007e-41 2.923504e-39    27
#> 2: HALLMARK_MYC_TARGETS_V1   35.385088 2.083184e-27 5.207959e-26    29
#> 3:    HALLMARK_E2F_TARGETS    4.219512 1.444580e-03 2.407633e-02     8
#>    gene_set_lengths target_set_lengths
#>               <num>              <int>
#> 1:               58                 50
#> 2:              200                 50
#> 3:              200                 50
```

To run the hypergeometric test over a list of target gene sets,
`bixverse` provides `gse_hypergeometric_list` that leverages
[rayon](https://github.com/rayon-rs/rayon) in Rust for multi-threading.
This makes this function very fast even on large gene set libraries with
lots of target gene sets.

``` r
target_genes_2 <- c(
  sample(h_gene_sets_ls[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]], 20),
  sample(h_gene_sets_ls[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]], 25)
)

target_list <- list(
  set_1 = target_genes_1,
  set_2 = target_genes_2,
  set_3 = c("random gene 1", "random gene 2", "random gene 3")
)

results_multiple <- gse_hypergeometric_list(
  target_genes_list = target_list,
  gene_set_list = h_gene_sets_ls
)

head(results_multiple)
#>    target_set_name odds_ratios        pvals          fdr  hits gene_set_lengths
#>             <char>       <num>        <num>        <num> <num>            <num>
#> 1:           set_2  162.370744 4.603263e-43 2.301631e-41    30               91
#> 2:           set_1  170.868035 5.740508e-41 2.870254e-39    27               58
#> 3:           set_1   35.410526 2.043327e-27 5.108317e-26    29              200
#> 4:           set_2   29.978469 5.916124e-22 1.479031e-20    24              200
#> 5:           set_2   12.001728 2.323959e-10 3.873265e-09    15              201
#> 6:           set_2    8.575841 2.503092e-07 3.128865e-06    12              200
#>                         gene_set_name target_set_lengths
#>                                <char>              <int>
#> 1:   HALLMARK_IL6_JAK_STAT3_SIGNALING                 45
#> 2:            HALLMARK_MYC_TARGETS_V2                 50
#> 3:            HALLMARK_MYC_TARGETS_V1                 50
#> 4:   HALLMARK_TNFA_SIGNALING_VIA_NFKB                 45
#> 5:     HALLMARK_INFLAMMATORY_RESPONSE                 45
#> 6: HALLMARK_INTERFERON_GAMMA_RESPONSE                 45
```

To exemplify the speed in Rust, in this case we run for 250 target gene
sets over-enrichment tests against 5000 gene sets sampled from 20k
genes. This totals 1.25m hypergeometric tests. Pending on your system,
you can expect this to be done in a few seconds.

``` r
seed = 10101L

set.seed(seed)

# Define the parameters
universe <- sprintf("gene_%i", 1:20000)
gene_sets_no <- 5000
target_gene_sets_no <- 250

# Generate random gene sets
gene_sets <- purrr::map(
  1:gene_sets_no,
  ~ {
    set.seed(seed + .x + 1)
    size <- sample(20:100, 1)
    sample(universe, size, replace = FALSE)
  }
)

names(gene_sets) <- purrr::map_chr(
  1:gene_sets_no,
  ~ {
    set.seed(seed + .x + 1)
    paste(sample(LETTERS, 3), collapse = "")
  }
)

# Generate random target sets
target_gene_sets <- purrr::map(
  1:target_gene_sets_no,
  ~ {
    set.seed(.x * seed)
    size <- sample(50:100, 1)
    sample(universe, size, replace = FALSE)
  }
)

names(target_gene_sets) <- purrr::map_chr(
  1:target_gene_sets_no,
  ~ {
    set.seed(seed + .x + 1)
    paste(sample(letters, 3), collapse = "")
  }
)

# Run the Rust function
tictoc::tic()
rs_results_example <- gse_hypergeometric_list(
  target_genes_list = target_gene_sets,
  gene_set_list = gene_sets
)
tictoc::toc()
#> 0.623 sec elapsed
```

### Gene ontology aware enrichment tests (for sets)

Additionally, `bixverse` provides Rust-accelerated functions to run an
ontology- aware enrichment test for the gene ontology in particular -
the concept is based on [Adrian et
al.](https://academic.oup.com/bioinformatics/article/22/13/1600/193669).
Briefly, the function will start in the leafs of the ontology (lowest
level) and calculate the enrichment statistic for that level. If an
enrichment test reaches a threshold, the genes from that term are
removed from all of its ancestor terms, ensuring that the more specific
terms have stronger test statistics. To make it easy for the users, we
provide the human gene ontology as data in the package. This is a
data.table with 5 columns:

- go_id: the gene ontology term identifier.
- go_name: the name gene ontology term name.
- ensembl_id: the ensembl identifiers belonging to this term.
- ancestors: in the ontology DAG all identified ancestors of this
  specific term
- depth: the depth of that term (identified via BFS) for this term.
  Higher values indicate further distance from the top terms.

These columns are expected by the `gene_ontology_data` class that stores
the data and exposes it in Rust.

``` r
# Load in the package provided data and process it into the right format
go_data_dt <- get_go_data_human()
#> Loading the data from the package.
#> Processing data for the gene_ontology class.
#> Warning: The `father` argument of `dfs()` is deprecated as of igraph 2.2.0.
#> ℹ Please use the `parent` argument instead.
#> ℹ The deprecated feature was likely used in the bixverse package.
#>   Please report the issue to the authors.

# Generate the S7 object for subsequent usage
go_data_s7 <- gene_ontology_data(go_data_dt, min_genes = 3L)
```

Similar to the more generic hypergeometric tests, you can provide either
a single set of a target genes or a list of target genes. Leveraging
Rusts borrowing, the overhead for multiple target gene sets is reduced
and hence very fast.

``` r
go_aware_res <- gse_go_elim_method(
  object = go_data_s7,
  target_genes = target_genes_1
)

head(go_aware_res)
#>         go_id                                     go_name odds_ratios
#>        <char>                                      <char>       <num>
#> 1: GO:0003723                                 RNA binding   14.812443
#> 2: GO:0032040                    small-subunit processome   44.202960
#> 3: GO:0042274          ribosomal small subunit biogenesis   43.541132
#> 4: GO:0005730                                   nucleolus    6.009426
#> 5: GO:0000055 ribosomal large subunit export from nucleus  455.891304
#> 6: GO:0031428    box C/D methylation guide snoRNP complex  455.891304
#>           pvals          fdr  hits gene_set_lengths
#>           <num>        <num> <num>            <num>
#> 1: 3.122680e-16 3.441817e-12    23             1205
#> 2: 1.623943e-08 6.489265e-05     6               72
#> 3: 1.766267e-08 6.489265e-05     6               73
#> 4: 8.523856e-08 2.348748e-04    18             1866
#> 5: 2.368116e-07 4.350230e-04     3                6
#> 6: 2.368116e-07 4.350230e-04     3                6
```

Here is the example on how to run this over a list.

``` r
go_aware_res_2 <- gse_go_elim_method_list(
  object = go_data_s7,
  target_gene_list = target_list
)

head(go_aware_res_2)
#>    target_set_name                                              go_name
#>             <char>                                               <char>
#> 1:           set_1                                          RNA binding
#> 2:           set_2 cell surface receptor signaling pathway via JAK-STAT
#> 3:           set_2                                inflammatory response
#> 4:           set_2                                    cytokine activity
#> 5:           set_2                     external side of plasma membrane
#> 6:           set_2                           cellular response to virus
#>         go_id odds_ratios        pvals          fdr  hits gene_set_lengths
#>        <char>       <num>        <num>        <num> <num>            <num>
#> 1: GO:0003723    14.81244 3.122680e-16 3.441817e-12    23             1205
#> 2: GO:0007259    97.16563 8.704402e-15 9.587028e-11     9               66
#> 3: GO:0006954    22.90818 1.408854e-13 7.758558e-10    14              447
#> 4: GO:0005125    31.85212 1.130962e-12 4.152138e-09    11              235
#> 5: GO:0009897    19.25374 1.941168e-10 5.133756e-07    11              379
#> 6: GO:0098586    56.46451 2.330559e-10 5.133756e-07     7               79
```

#### Speed of the functions

Also this functions were designed with speed in mind (leveraging Rust’s
borrowing and parallelisations where possible). To run this over 100
potential target gene sets takes only a few seconds.

``` r
go_gene_universe <- unique(unlist(go_data_dt$ensembl_id))

# Generate random target sets

go_target_sets_no <- 100L

seed <- 246L

go_target_gene_sets <- purrr::map(
  1:go_target_sets_no,
  ~ {
    set.seed(.x * seed)
    size <- sample(50:100, 1)
    sample(go_gene_universe, size, replace = FALSE)
  }
)

names(go_target_gene_sets) <- purrr::map_chr(
  1:go_target_sets_no,
  ~ {
    set.seed(seed + .x + 1)
    paste(sample(letters, 3), collapse = "")
  }
)

# Quick test with tictoc
tictoc::tic()
rs_results_example <- gse_go_elim_method_list(
  object = go_data_s7,
  target_gene_list = go_target_gene_sets
)
tictoc::toc()
#> 1.477 sec elapsed
```

### Alternative: simplifying results

Another approach is to run the enrichment over the gene ontology as is
and simplify subsequently. `bixverse` also has options for that. Let’s
run an example with the gene ontology:

``` r
# load the go data that is part of the package
go_data <- load_go_human_data()

# prepare the go <> gene information
min_genes <- 3L

go_genes <- go_data$go_to_genes
go_genes_ls <- split(go_genes$ensembl_id, go_genes$go_id)
go_genes_ls <- purrr::keep(go_genes_ls, \(x) length(x) > min_genes)

# run the hypergeometric test
go_results_unfiltered <- gse_hypergeometric(
  target_genes = target_genes_1,
  gene_set_list = go_genes_ls
)

head(go_results_unfiltered)
#>    gene_set_name odds_ratios        pvals          fdr  hits gene_set_lengths
#>           <char>       <num>        <num>        <num> <num>            <num>
#> 1:    GO:0003723   22.120609 1.723343e-23 1.579099e-19    31             1546
#> 2:    GO:0006364   35.837205 2.936032e-12 1.345143e-08    10              159
#> 3:    GO:0005730    8.860092 6.601460e-12 1.581347e-08    23             1927
#> 4:    GO:0042254   43.069954 6.903186e-12 1.581347e-08     9              118
#> 5:    GO:0005634    8.314215 9.730968e-12 1.783297e-08    39             6736
#> 6:    GO:0005654    6.752334 9.279447e-11 1.417126e-07    30             4005
#>    target_set_lengths
#>                 <int>
#> 1:                 50
#> 2:                 50
#> 3:                 50
#> 4:                 50
#> 5:                 50
#> 6:                 50
```

Now that we have the results, we can simplify them. The approach here is
to identify similar terms via the Wang similarity and within a subset of
similar terms take the term with the best test statistic, i.e., lowest
FDR. Should the terms have the same FDR, the function will select for
the most specific term (i.e., lower in the ontology).

``` r
go_parent_child_dt <- go_data$gene_ontology[
  relationship %in% c("is_a", "part_of")
] %>%
  setnames(
    old = c("from", "to", "relationship"),
    new = c("parent", "child", "type")
  )

go_results_simplified <- simplify_hypergeom_res(
  res = go_results_unfiltered,
  parent_child_dt = go_parent_child_dt,
  weights = setNames(c(0.8, 0.6), c("is_a", "part_of"))
)

head(go_results_simplified)
#>    gene_set_name odds_ratios        pvals          fdr  hits gene_set_lengths
#>           <char>       <num>        <num>        <num> <num>            <num>
#> 1:    GO:0003723   22.120609 1.723343e-23 1.579099e-19    31             1546
#> 2:    GO:0006364   35.837205 2.936032e-12 1.345143e-08    10              159
#> 3:    GO:0005730    8.860092 6.601460e-12 1.581347e-08    23             1927
#> 4:    GO:0042254   43.069954 6.903186e-12 1.581347e-08     9              118
#> 5:    GO:0032040   44.202960 1.623943e-08 2.023038e-05     6               72
#> 6:    GO:0000055  455.891304 2.368116e-07 2.169905e-04     3                6
#>    target_set_lengths
#>                 <int>
#> 1:                 50
#> 2:                 50
#> 3:                 50
#> 4:                 50
#> 5:                 50
#> 6:                 50
```

### GSEA

`bixverse` also has implementations of the GSEA based on [Subramanian et
al](https://www.pnas.org/doi/10.1073/pnas.0506580102) to test against
continuous scores. The package provides also the fgsea implementations
of [Korotkevich et
al.](https://www.biorxiv.org/content/10.1101/060012v3). These can be run
like:

``` r
# This is the example of fgsea

library("fgsea")

data(examplePathways)
data(exampleRanks)

set.seed(42L)

fgsea_res <- fgsea(
  pathways = examplePathways, 
  stats = exampleRanks,
  minSize = 15,
  maxSize = 500
) %>% setorder(pathway)

head(fgsea_res)
#>                                                                                    pathway
#>                                                                                     <char>
#> 1:                                                                1221633_Meiotic_Synapsis
#> 2:                                   1445146_Translocation_of_Glut4_to_the_Plasma_Membrane
#> 3: 442533_Transcriptional_Regulation_of_Adipocyte_Differentiation_in_3T3-L1_Pre-adipocytes
#> 4:                                                                  508751_Circadian_Clock
#> 5:                                               5334727_Mus_musculus_biological_processes
#> 6:                                        573389_NoRC_negatively_regulates_rRNA_expression
#>         pval      padj    log2err         ES        NES  size
#>        <num>     <num>      <num>      <num>      <num> <int>
#> 1: 0.5490534 0.7262873 0.06674261  0.2885754  0.9399884    27
#> 2: 0.6952862 0.8366277 0.05445560  0.2387284  0.8366856    39
#> 3: 0.1122449 0.2599823 0.21392786 -0.3640706 -1.3460572    31
#> 4: 0.7826888 0.8799951 0.05312981  0.2516324  0.7287088    17
#> 5: 0.3580060 0.5579562 0.08197788  0.2469065  1.0498921   106
#> 6: 0.4198895 0.6197865 0.08407456  0.3607407  1.0446784    17
#>                                 leadingEdge
#>                                      <list>
#> 1:                        15270,12189,71846
#> 2:  17918,19341,20336,22628,22627,20619,...
#> 3: 76199,19014,26896,229003,17977,17978,...
#> 4:                              20893,59027
#> 5:  60406,19361,15270,20893,12189,68240,...
#> 6:                              60406,20018
```

And this is how you run it via `bixverse`.

``` r
bixverse_fgsea <- calc_fgsea(
  stats = exampleRanks,
  pathways = examplePathways,
  gsea_params = params_gsea(min_size = 15L)
) %>% setorder(pathway_name)

head(bixverse_fgsea)
#>            es        nes     pvals n_more_extreme  size
#>         <num>      <num>     <num>          <num> <num>
#> 1:  0.2885754  0.9321129 0.5588723            336    27
#> 2:  0.2387284  0.8353752 0.7140600            451    39
#> 3: -0.3640706 -1.3241710 0.1282051             49    31
#> 4:  0.2516324  0.7211793 0.8109541            458    17
#> 5:  0.2469065  1.0551934 0.3751804            259   106
#> 6:  0.3607407  1.0338840 0.4204947            237    17
#>                                                                               pathway_name
#>                                                                                     <char>
#> 1:                                                                1221633_Meiotic_Synapsis
#> 2:                                   1445146_Translocation_of_Glut4_to_the_Plasma_Membrane
#> 3: 442533_Transcriptional_Regulation_of_Adipocyte_Differentiation_in_3T3-L1_Pre-adipocytes
#> 4:                                                                  508751_Circadian_Clock
#> 5:                                               5334727_Mus_musculus_biological_processes
#> 6:                                        573389_NoRC_negatively_regulates_rRNA_expression
#>                                leading_edge    log2err       fdr
#>                                      <list>      <num>     <num>
#> 1:                  15270,12189,71846,19357 0.06407038 0.7359532
#> 2:  17918,19341,20336,22628,22627,20619,... 0.05029481 0.8557113
#> 3: 76199,19014,26896,229003,17977,17978,... 0.19991523 0.2900703
#> 4:                        20893,59027,19883 0.04959020 0.9063587
#> 5:  60406,19361,15270,20893,12189,68240,... 0.07707367 0.5816288
#> 6:                 60406,20018,245688,20017 0.08175156 0.6207015
```

We can appreciate the (more or less) same p-values estimated by both
methods.

``` r
plot(
  x = -log10(fgsea_res$pval),
  y = -log10(bixverse_fgsea$pvals),
  xlab = "-log10(pval) fgsea",
  ylab = "-log10(pval) bixverse",
  main = "fgsea and bixverse"
)
```

![](gse_methods_files/figure-html/p-value%20comparison%20between%20fgsea%20and%20bixverse%20fgsea-1.png)

The speed is very comparable with the Rcpp implementation of fgsea,
indicating no issues in terms of performance.

``` r
microbenchmark::microbenchmark(
  fgsea = fgsea(
    pathways = examplePathways, 
    stats = exampleRanks,
    minSize = 15,
    maxSize = 500
  ),
  rust = calc_fgsea(
    stats = exampleRanks,
    pathways = examplePathways,
    gsea_params = params_gsea(min_size = 15L)
  ),
  times = 10L
)
#> Unit: seconds
#>   expr      min       lq     mean   median       uq      max neval
#>  fgsea 2.262689 2.361076 2.573897 2.409774 2.849547 3.072855    10
#>   rust 2.496505 2.500094 2.513522 2.509437 2.522316 2.548479    10
```

### GSEA gene ontology aware

`bixverse` also implements GSEA in a gene-ontology aware form. In this
case, you use the same object as for the hypergeometric test, but
instead provide the ranked list. Similar to the hyper geometric version,
the ontology is traversed from the lowest level first and the
permutation-based p-values are calculated for a given term. Should this
p-value be below the elimination threshold, the genes from that given
gene ontology term are removed from all of its ancestors, yielding more
specific GO terms. Once this has been done with the full ontology, the
multi-level method from fgsea is used to estimate lower p-values than
just possible via permutation for terms that show significance.

``` r
gene_universe_go <- unique(
  unlist(go_data_dt[, "ensembl_id"], use.names = FALSE)
)

set.seed(42L)

random_stats <- rnorm(length(gene_universe_go))
names(random_stats) <- gene_universe_go

go_gsea_res <- fgsea_go_elim(
  object = go_data_s7,
  stats = random_stats
)

head(go_gsea_res)
#>         go_id         es       nes  size        pvals n_more_extreme
#>        <char>      <num>     <num> <num>        <num>          <num>
#> 1: GO:0031465 -0.9025468 -1.925965     6 6.165332e-05              0
#> 2: GO:0050877  0.5116421  1.950886    44 2.410338e-04              2
#> 3: GO:0072675 -0.8744037 -1.865910     6 2.625170e-04              0
#> 4: GO:0006066 -0.6718553 -2.026201    18 2.712234e-04              0
#> 5: GO:0071392 -0.5217830 -1.904151    36 3.788231e-04              0
#> 6: GO:0002090 -0.8302363 -1.868975     7 6.827804e-04              3
#>                                                                                           leading_edge
#>                                                                                                 <list>
#> 1:                                                     ENSG00000110844,ENSG00000180370,ENSG00000160685
#> 2: ENSG00000131263,ENSG00000229674,ENSG00000172987,ENSG00000108878,ENSG00000050165,ENSG00000143061,...
#> 3:                                                                                     ENSG00000170613
#> 4: ENSG00000287395,ENSG00000139354,ENSG00000211659,ENSG00000160345,ENSG00000007237,ENSG00000139679,...
#> 5: ENSG00000118193,ENSG00000134531,ENSG00000166897,ENSG00000162980,ENSG00000170348,ENSG00000103018,...
#> 6:                                     ENSG00000061938,ENSG00000169507,ENSG00000289721,ENSG00000083814
#>      log2err       fdr                                 go_name
#>        <num>     <num>                                  <char>
#> 1: 0.5384341 0.4815124  Cul4B-RING E3 ubiquitin ligase complex
#> 2: 0.5188481 0.5295637                  nervous system process
#> 3: 0.4984931 0.5295637                       osteoclast fusion
#> 4: 0.4984931 0.5295637               alcohol metabolic process
#> 5: 0.4984931 0.5917217 cellular response to estradiol stimulus
#> 6: 0.4772708 0.8290480  regulation of receptor internalization
```
