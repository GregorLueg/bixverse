# Thinking about single cells the bixverse way

## Intro

Before diving into this vignette, it would be really recommended to read
the [design
choices](https://gregorlueg.github.io/bixverse/articles/design_single_cell.html)
to understand **why** single cell has been implemented the way it is in
this package. This will make understanding the downstream text much
easier. The aim of this vignette is to take some small synthetic data
and point out some quirks due to design choices made.

``` r
library(bixverse)
```

### Synthetic data

Let’s start with some synthetic data that is powering all the tests in
the package. It’s very simple data with 1,000 cells, 100 genes and three
very cleary defined cell groups. Genes 1 to 10 defined cell_type 1,
genes 11 to 20 the next cell type and 21 to 30 the last cell type. The
remaining genes are noise genes.

``` r
single_cell_test_data <- generate_single_cell_test_data()

str(single_cell_test_data)
#> List of 3
#>  $ counts:Formal class 'dgRMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ p       : int [1:1001] 0 61 115 174 233 295 359 420 484 542 ...
#>   .. ..@ j       : int [1:60179] 0 1 2 3 4 5 6 7 8 9 ...
#>   .. ..@ Dim     : int [1:2] 1000 100
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : chr [1:1000] "cell_0001" "cell_0002" "cell_0003" "cell_0004" ...
#>   .. .. ..$ : chr [1:100] "gene_001" "gene_002" "gene_003" "gene_004" ...
#>   .. ..@ x       : num [1:60179] 14 24 7 13 9 44 23 46 6 1 ...
#>   .. ..@ factors : list()
#>  $ obs   :Classes 'data.table' and 'data.frame': 1000 obs. of  3 variables:
#>   ..$ cell_id    : chr [1:1000] "cell_0001" "cell_0002" "cell_0003" "cell_0004" ...
#>   ..$ cell_grp   : chr [1:1000] "cell_type_1" "cell_type_2" "cell_type_3" "cell_type_1" ...
#>   ..$ batch_index: num [1:1000] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ var   :Classes 'data.table' and 'data.frame': 100 obs. of  2 variables:
#>   ..$ gene_id   : chr [1:100] "gene_001" "gene_002" "gene_003" "gene_004" ...
#>   ..$ ensembl_id: chr [1:100] "ens_001" "ens_002" "ens_003" "ens_004" ...
#>   ..- attr(*, ".internal.selfref")=<externalptr>
```

We have a count matrix with pseudo raw counts, an obs table and a var
table. Let us load this into a `bixverse` SingleCell class and explore
what’s actually happening.

### Loading in data

Initialising an experiment works like this. We will not save the data
into any persistent storage, but just leave it as is… If you read the
design docs, you know why you have to provide a path here. If you do not
know why, go back to the docs.

``` r
test_dir <- file.path(tempdir(), "single_cell_test")
dir.create(test_dir, showWarnings = FALSE, recursive = TRUE)

sc_object <- SingleCells(dir_data = test_dir)

sc_object
#> Single cell experiment (Single Cells).
#>   No cells (original): 0
#>    To keep n: 0
#>   No genes: 0
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

Is something in the directory?

``` r
list.files(test_dir)
#> character(0)
```

Nope. But let’s now load in the data and see what happens:

``` r
sc_object <- load_r_data(
  object = sc_object,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  var = single_cell_test_data$var,
  # we let all cells pass for now
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 0L,
    min_lib_size = 0L,
    min_cells = 0L,
    target_size = 1e3 # target size of 1000
  ),
  streaming = FALSE,
  .verbose = TRUE
)
#> Writing counts to disk.
#> Generating gene-based data.
#> Writing to the DuckDB.
#> Setting internal mapping.

list.files(test_dir)
#> [1] "counts_cells.bin" "counts_genes.bin" "sc_duckdb.db"
```

So, what has just happened … ? Actually quite a lot of things:

1.  The R matrix is being pulled into Rust and checked against the
    number of min_unique_genes, the min_lib_size and in how many cells a
    given gene needs to expressed to be included. Once the genes and
    cells have been identified that pass the minimum threshold (in that
    order: first gene filtering -\> then check the cells against that
    reduced feature space), the cells are written to the first binary
    file: `counts_cells.bin`. At the same time the counts are
    normalised: division by library size, projected against desired
    target size and log transformed.
2.  Rust now transposes under the hood the data that is at the moment
    very nice to read in from a cell-centric level, less so from a
    gene-centric level. That is why we also have `counts_genes.bin`.
3.  Lastly, the observation and variable data is being reduced based on
    step 1 and written to the `sc_duckdb.db` and stored there.

If we now print the class:

``` r
sc_object
#> Single cell experiment (Single Cells).
#>   No cells (original): 1000
#>    To keep n: 1000
#>   No genes: 100
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

We can now see the number of cells and genes stored in the object.

### Interacting with the data

Let’s now look into how we interact with the data… We will actually
apply some filters here (the synthetic data is weird, but does the job
of assumption testing).

``` r
sc_object <- load_r_data(
  object = sc_object,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  var = single_cell_test_data$var,
  # we let all cells pass for now
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 45L,
    min_lib_size = 300L,
    min_cells = 500L, # the data is weird and a lot cells express genes
    target_size = 1e3 # target size of 1000
  ),
  streaming = FALSE,
  .verbose = TRUE
)
#> Writing counts to disk.
#> Generating gene-based data.
#> Writing to the DuckDB.
#> Setting internal mapping.

sc_object
#> Single cell experiment (Single Cells).
#>   No cells (original): 929
#>    To keep n: 929
#>   No genes: 81
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

We now have overwritten the files on disk for the filtered versions. 929
cells and 81 genes pass the threshold. But how do I interact with the
data?

#### Understanding mapping and filtering

As data is stored on disk, we need to carefully synchronise state here
between the DuckDB, the binary files and tell the Rust kernel which
cells to use and what to do… There is a bunch of helpers for this. You
might have realised in the print like `To keep n`. This is very
important! Let’s check this out

``` r
get_cells_to_keep(sc_object)[1:10]
#>  [1] 0 1 2 3 4 5 6 7 8 9
```

What is cells to keep? Well, these are the (0-based) indices of cells
that shall be included in the analysis, think HVG selection, PCA, kNN
generation, etc. This is **THE** foot gun of this package. The
[`get_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/get_cells_to_keep.md)
indices will be the ones loaded in by Rust. The other cells are not
lost, they are still on disk (a massive difference between filtering in
memory and on-disk work!). Also the DuckDB state will get updated, but
let’s get to this later. We also have other getters.

``` r
# get the cell identifiers - filtered for genes to keep; if set to FALSE
# everything will be returned
print(get_cell_names(sc_object, filtered = TRUE)[1:10])
#>  [1] "cell_0001" "cell_0003" "cell_0004" "cell_0005" "cell_0006" "cell_0007"
#>  [7] "cell_0008" "cell_0009" "cell_0010" "cell_0011"

print(get_gene_names(sc_object)[1:10])
#>  [1] "gene_001" "gene_002" "gene_003" "gene_004" "gene_005" "gene_006"
#>  [7] "gene_007" "gene_008" "gene_009" "gene_010"
```

Let’s quickly check out how you access the obs table and then let’s
exemplify THE footgun of the paper:

``` r
# a more R-native way
sc_object[[]]
#>      cell_idx   cell_id    cell_grp batch_index   nnz lib_size to_keep
#>         <int>    <char>      <char>       <num> <num>    <num>  <lgcl>
#>   1:        1 cell_0001 cell_type_1           1    52      371    TRUE
#>   2:        2 cell_0003 cell_type_3           1    52      495    TRUE
#>   3:        3 cell_0004 cell_type_1           1    48      551    TRUE
#>   4:        4 cell_0005 cell_type_2           1    53      422    TRUE
#>   5:        5 cell_0006 cell_type_3           1    57      394    TRUE
#>  ---                                                                  
#> 925:      925 cell_0995 cell_type_2           1    53      548    TRUE
#> 926:      926 cell_0996 cell_type_3           1    57      545    TRUE
#> 927:      927 cell_0997 cell_type_1           1    56      543    TRUE
#> 928:      928 cell_0998 cell_type_2           1    49      487    TRUE
#> 929:      929 cell_1000 cell_type_1           1    60      474    TRUE
```

``` r
# a more R-native way: only first three rows
sc_object[[1:3L]]
#>    cell_idx   cell_id    cell_grp batch_index   nnz lib_size to_keep
#>       <int>    <char>      <char>       <num> <num>    <num>  <lgcl>
#> 1:        1 cell_0001 cell_type_1           1    52      371    TRUE
#> 2:        2 cell_0003 cell_type_3           1    52      495    TRUE
#> 3:        3 cell_0004 cell_type_1           1    48      551    TRUE
```

``` r
# a more R-native way: specific columns
sc_object[[c("cell_idx", "cell_id", "cell_grp")]]
#>      cell_idx   cell_id    cell_grp
#>         <int>    <char>      <char>
#>   1:        1 cell_0001 cell_type_1
#>   2:        2 cell_0003 cell_type_3
#>   3:        3 cell_0004 cell_type_1
#>   4:        4 cell_0005 cell_type_2
#>   5:        5 cell_0006 cell_type_3
#>  ---                               
#> 925:      925 cell_0995 cell_type_2
#> 926:      926 cell_0996 cell_type_3
#> 927:      927 cell_0997 cell_type_1
#> 928:      928 cell_0998 cell_type_2
#> 929:      929 cell_1000 cell_type_1
```

All of these ways above call R primitives. There is also a more general
getter here:

``` r
get_sc_obs(sc_object, filtered = FALSE)
#>      cell_idx   cell_id    cell_grp batch_index   nnz lib_size to_keep
#>         <int>    <char>      <char>       <num> <num>    <num>  <lgcl>
#>   1:        1 cell_0001 cell_type_1           1    52      371    TRUE
#>   2:        2 cell_0003 cell_type_3           1    52      495    TRUE
#>   3:        3 cell_0004 cell_type_1           1    48      551    TRUE
#>   4:        4 cell_0005 cell_type_2           1    53      422    TRUE
#>   5:        5 cell_0006 cell_type_3           1    57      394    TRUE
#>  ---                                                                  
#> 925:      925 cell_0995 cell_type_2           1    53      548    TRUE
#> 926:      926 cell_0996 cell_type_3           1    57      545    TRUE
#> 927:      927 cell_0997 cell_type_1           1    56      543    TRUE
#> 928:      928 cell_0998 cell_type_2           1    49      487    TRUE
#> 929:      929 cell_1000 cell_type_1           1    60      474    TRUE
# this one can retrieve unfiltered cells!
```

Let’s pretend we want to do mitochondrial (or other gene set-based
removal)

``` r
# pseudo gene sets
gs_of_interest <- list(
  gs_1 = c("gene_001", "gene_002", "gene_003", "gene_004"),
  gs_2 = c("gene_096", "gene_097", "gene_100")
)

# calculate the gene set proportions
sc_object <- gene_set_proportions_sc(
  sc_object,
  gs_of_interest,
  .verbose = TRUE
)
```

If we now look at the obs again:

``` r
sc_object[[]]
#>      cell_idx   cell_id    cell_grp batch_index   nnz lib_size to_keep
#>         <int>    <char>      <char>       <num> <num>    <num>  <lgcl>
#>   1:        1 cell_0001 cell_type_1           1    52      371    TRUE
#>   2:        2 cell_0003 cell_type_3           1    52      495    TRUE
#>   3:        3 cell_0004 cell_type_1           1    48      551    TRUE
#>   4:        4 cell_0005 cell_type_2           1    53      422    TRUE
#>   5:        5 cell_0006 cell_type_3           1    57      394    TRUE
#>  ---                                                                  
#> 925:      925 cell_0995 cell_type_2           1    53      548    TRUE
#> 926:      926 cell_0996 cell_type_3           1    57      545    TRUE
#> 927:      927 cell_0997 cell_type_1           1    56      543    TRUE
#> 928:      928 cell_0998 cell_type_2           1    49      487    TRUE
#> 929:      929 cell_1000 cell_type_1           1    60      474    TRUE
#>             gs_1        gs_2
#>            <num>       <num>
#>   1: 0.156334236 0.008086253
#>   2: 0.022222223 0.006060606
#>   3: 0.181488201 0.000000000
#>   4: 0.007109005 0.033175357
#>   5: 0.017766498 0.015228426
#>  ---                        
#> 925: 0.021897810 0.041970804
#> 926: 0.007339450 0.034862384
#> 927: 0.267034978 0.173112333
#> 928: 0.016427105 0.047227927
#> 929: 0.187763706 0.040084388
```

We now have two new columns here! That worked well… Let’s now do
filtering to really understand what’s happening under the hood:

``` r
threshold <- 0.05

cells_to_keep <- sc_object[[]][gs_2 < threshold, cell_id]

sc_object <- set_cells_to_keep(sc_object, cells_to_keep)

sc_object
#> Single cell experiment (Single Cells).
#>   No cells (original): 929
#>    To keep n: 682
#>   No genes: 81
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

The print changed now! We have now told the class that only a subset of
682 cells shall be used in any downstream analyses. This is very
important! Most functions will only return and work on the filtered
cells, and you have to carefully think and remember what you did when.
Also:

``` r
obs_unfiltered <- get_sc_obs(sc_object, filtered = FALSE)

obs_filtered <- get_sc_obs(sc_object, filtered = TRUE)

obs_r_primitive <- sc_object[[]]

cat(sprintf("No of rows unfiltered: %i\n", nrow(obs_unfiltered)))
#> No of rows unfiltered: 929
cat(sprintf("No of rows filtered: %i\n", nrow(obs_filtered)))
#> No of rows filtered: 682
cat(sprintf("No of rows from R primitive: %i\n", nrow(obs_r_primitive)))
#> No of rows from R primitive: 682
```

If you understand the concept shown above, you are good to go to analyse
millions of cells on small compute without having 100’s of GB of memory
available. Let’s just talk also about the 2nd big foot gun: **the
counts**.

#### Understanding count retrieval

You have actually two files on disk… If you read the design doc you know
why (if you don’t, please read it!). Now let’s look into what it
ACTUALLY means… Let’s say you are interested in the first 50 cells.

``` r
raw_counts_from_cells <- sc_object[
  1:50L,
  ,
  assay = "raw",
  return_format = "cell"
]

print(class(raw_counts_from_cells))
#> [1] "dgRMatrix"
#> attr(,"package")
#> [1] "Matrix"

raw_counts_from_genes <- sc_object[
  1:50L,
  ,
  assay = "raw",
  return_format = "gene"
]

print(class(raw_counts_from_genes))
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"
```

What is happening here? In both cases, we are pulling the raw counts. In
the first case, we do it smartly… Rows -\> cells, so, we are telling R
to tell Rust to retrieve this as CSR. Which means we ONLY load the cells
in that are needed. The second approach is the stupid way… We are
telling R to tell Rust to load in ALL genes and then filter to the
desired cells. The impact of this can be seen here:

``` r
microbenchmark::microbenchmark(
  the_correct_way = {
    sc_object[
      1:50L,
      ,
      assay = "raw",
      return_format = "cell"
    ]
  },
  the_incorrect_way = {
    sc_object[
      1:50L,
      ,
      assay = "raw",
      return_format = "gene"
    ]
  },
  times = 10L
)
#> Unit: milliseconds
#>               expr      min       lq     mean   median       uq      max neval
#>    the_correct_way 1.273197 1.281084 1.363823 1.350251 1.430671 1.539825    10
#>  the_incorrect_way 2.020301 2.148822 2.289115 2.238745 2.323718 3.027833    10
```

The difference seems marginal here, but it WILL bite you if you do this
on large data sets. Also, the memory advantages this approach gives you
usually will disappear if you force the wrong format. Also, another
thing to REALLY consider… The indices you supply in R should match what
the Rust files have and again you need to THINK about the cells to keep.
A nasty example here:

``` r
counts_filtered <- get_sc_counts(
  object = sc_object,
  assay = "raw",
  return_format = "gene",
  # get from the R perspective cells 1:10 that were written to disk
  cell_indices = 1:10L,
  # get from the R perspective genes 1:10 that were written to disk
  gene_indices = 1:10L,
  use_cells_to_keep = TRUE
)

counts_unfiltered <- get_sc_counts(
  object = sc_object,
  assay = "raw",
  return_format = "gene",
  # get from the R perspective cells 1:10 that were written to disk
  cell_indices = 1:10L,
  # get from the R perspective genes 1:10 that were written to disk
  gene_indices = 1:10L,
  use_cells_to_keep = FALSE
)

cat(sprintf(
  "Number of cells: %i | Number of genes: %i (filtered version)\n",
  nrow(counts_filtered),
  ncol(counts_filtered)
))
#> Number of cells: 8 | Number of genes: 10 (filtered version)
cat(sprintf(
  "Number of cells: %i | Number of genes: %i (unfiltered version)\n",
  nrow(counts_unfiltered),
  ncol(counts_unfiltered)
))
#> Number of cells: 10 | Number of genes: 10 (unfiltered version)
```

Wait, what happened here … ? Because we set the cells to keep in the
step with the pseudo mitochondrial counts, these indices are now “bad
cells”. **They are not filtered, they are just not returned.** You can
control this also via the primitive syntax.

``` r
filtered_counts <- sc_object[1:10L, 1:10L, use_cells_to_keep = TRUE]

unfiltered_counts <- sc_object[1:10L, 1:10L, use_cells_to_keep = FALSE]

cat(sprintf(
  "Number of cells: %i | Number of genes: %i (filtered version - primitive)\n",
  nrow(counts_filtered),
  ncol(counts_filtered)
))
#> Number of cells: 8 | Number of genes: 10 (filtered version - primitive)
cat(sprintf(
  "Number of cells: %i | Number of genes: %i (unfiltered version - primitive)\n",
  nrow(counts_unfiltered),
  ncol(counts_unfiltered)
))
#> Number of cells: 10 | Number of genes: 10 (unfiltered version - primitive)
```

If you want to avoid thinking about indices altogether (which is
recommended), just use:

``` r
# just use cell labels and genes directly...
sc_object[
  c("cell_0001", "cell_0003", "cell_0100", "cell_0117"),
  c("gene_001", "gene_002", "gene_003", "gene_081")
]
#> 4 x 4 sparse Matrix of class "dgRMatrix"
#>           gene_001 gene_002 gene_003 gene_081
#> cell_0001       14       24        7        3
#> cell_0003        3        3        1        .
#> cell_0100        3       40       18        9
#> cell_0117        5        3        1        9
```

And that’s it. If you understand this vignette and the design choices,
we have answered the *why*, the *how* and can move on to the *what* in
the other vignettes with some example analyses and how to actually use
the package to analyse single cell experiments.
