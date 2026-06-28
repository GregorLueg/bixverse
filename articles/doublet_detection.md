# Doublet detection with bixverse

## Intro

Doublets - droplets containing two or more cells - are a common artefact
in droplet-based single cell experiments. If left undetected, they
masquerade as intermediate cell states or novel populations and quietly
corrupt downstream analyses. This vignette demonstrates three doublet
detection approaches available in `bixverse` and benchmarks them against
ground truth calls from [demuxlet](https://github.com/statgen/demuxlet),
which identifies doublets via genetic variation across pooled donors.

`bixverse` provides three computational doublet detection methods:

- **Scrublet** is a simulation-based approach ([Wolock, et al.,
  2019](https://doi.org/10.1016/j.cels.2018.11.005)). It generates
  synthetic doublets by averaging pairs of randomly selected cells,
  embeds both real and synthetic cells together, and scores each real
  cell by the proportion of its nearest neighbours that are synthetic. A
  bimodal score distribution is expected: singlets cluster near zero,
  doublets near one. A threshold on this score determines the final
  calls.

- **Boosted doublet detection** takes a different angle. Rather than
  simulating doublets, it trains a gradient-boosted tree ensemble to
  distinguish real cells from synthetic doublets generated in the same
  way as Scrublet. The classifier produces a probability per cell, and
  doublet calls are derived from that probability. This can be more
  robust when the score distribution from Scrublet is not cleanly
  bimodal, though it introduces its own hyperparameters (and is
  substantially slower).

- **scDblFinder** is a port of the
  [scDblFinder](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html)
  method from [Germain et
  al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9204188/), with one
  modification: it uses lightGBM with the histogram trick on the feature
  space instead of xgboost, which tends to be faster to fit. Like the
  boosted approach it trains a classifier on real versus synthetic
  doublets, but it combines the classifier output with cxds
  co-expression scores ([Bais and
  Kostka](https://pmc.ncbi.nlm.nih.gov/articles/PMC7703774/)) into a
  weighted final score. Across published benchmarks it tends to be the
  strongest of the three.

All three methods are purely computational and work without genetic or
hashing information. They complement experimental approaches like
demuxlet, which require pooled donors or hashtag antibodies. They also
have a failure mode if a large number of doublets dominates the original
data.

``` r

library(bixverse)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
```

## Loading the data

We use a PBMC data set with demuxlet ground truth calls. The demuxlet
output classifies each barcode as a singlet (`SNG`) or doublet (`DBL`)
(or ambiguous, but irrelevant for this situation), giving us something
to benchmark against.

``` r

doublet_path <- download_demuxlet_pbmc()

tempdir_doublet <- file.path(tempdir(), "demuxlet_bixverse")
dir.create(tempdir_doublet, showWarnings = FALSE, recursive = TRUE)

demuxlet_data <- fread(file.path(doublet_path, "demuxlet_calls.tsv"))

table(demuxlet_data$Call)
#> 
#>   AMB   DBL   SNG 
#>    24  1565 13030
```

Load in the data into bixverse

``` r

sc_object <- SingleCells(dir_data = tempdir_doublet)

mtx_io_params <- get_cell_ranger_params(doublet_path)
mtx_io_params$cells_as_rows <- TRUE

sc_object <- load_mtx(
  object = sc_object,
  sc_mtx_io_param = mtx_io_params,
  streaming = 0L
)
#>  Loading data directly into memory for CSR to CSC conversion.
#> Loading observations data from flat file into the DuckDB.
#> Loading variable data from flat file into the DuckDB.

sc_object
#> Single cell experiment (Single Cells).
#>   No cells (original): 14528
#>    To keep n: 14528
#>   No genes: 12622
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

We also define a small helper to compute precision, recall and F1
against the demuxlet calls.

``` r

doublet_metrics <- function(
  predicted,
  actual,
  pos_predicted = TRUE,
  pos_actual = "DBL"
) {
  tp <- sum(predicted == pos_predicted & actual == pos_actual)
  fp <- sum(predicted == pos_predicted & actual != pos_actual)
  fn <- sum(predicted != pos_predicted & actual == pos_actual)
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1 <- 2 * (precision * recall) / (precision + recall)
  list(precision = precision, recall = recall, f1 = f1)
}
```

## Doublet detection methods

Let’s check the two methods out and see how they behave

### Scrublet

``` r

scrublet <- scrublet_sc(
  object = sc_object,
  scrublet_params = params_scrublet(expected_doublet_rate = 0.12)
)

scrublet
#> ScrubletRes: 14528 cells, 1802 doublets (12.4%)
#>   Threshold:              0.2051
#>   Detected doublet rate:  12.4%
#>   Detectable fraction:    88.2%
#>   Overall doublet rate:   14.1%
#>   Simulated doublets:     21792
```

The `plot` method shows the score distribution for observed and
simulated cells. A clean separation between the two modes indicates the
threshold is sensible.

``` r

plot(scrublet)
```

![](doublet_detection_files/figure-html/scrublet%20plot-1.png)

How does the automatic thresholding perform against the demuxlet ground
truth?

``` r

scrublet_dt <- get_data(scrublet)
scrublet_dt[, Barcode := get_cell_names(sc_object)]
scrublet_dt <- merge(scrublet_dt, demuxlet_data, by = "Barcode")

doublet_metrics(predicted = scrublet_dt$doublet, actual = scrublet_dt$Call)
#> $precision
#> [1] 0.6553829
#> 
#> $recall
#> [1] 0.7575369
#> 
#> $f1
#> [1] 0.702767
```

If you want to add the information to the DuckDB within the object, you
can just run:

``` r

sc_object <- add_sc_new_obs(
  object = sc_object,
  obs_data = get_data(scrublet)
)

head(sc_object)
#>    cell_idx          cell_id   nnz lib_size to_keep doublet doublet_score
#>       <num>           <char> <num>    <num>  <lgcl>  <lgcl>         <num>
#> 1:        1 AAACATACAATGCC-1   851     1313    TRUE    TRUE  0.3265306652
#> 2:        2 AAACATACATTTCC-1   876     1537    TRUE   FALSE  0.0668662637
#> 3:        3 AAACATACCAGAAA-1   713     1149    TRUE   FALSE  0.0209840816
#> 4:        4 AAACATACCAGCTA-1   948     1588    TRUE    TRUE  0.3494623899
#> 5:        5 AAACATACCATGCA-1   337      499    TRUE   FALSE  0.0006016847
#> 6:        6 AAACATACCTCGCT-1   849     1494    TRUE   FALSE  0.0685483888
```

As we can see with added now the Scrublet data

#### Manual threshold adjustment

The automatic threshold does not always land in the right place. If the
score distribution suggests a different cut-off, `call_doublets_manual`
lets you override it.

``` r

scrublet_adj <- call_doublets_manual(scrublet, threshold = 0.25)
#> Detected doublet rate = 11.3%
#> Estimated detectable doublet fraction = 84.3%
#> Overall doublet rate:
#>   Estimated = 13.4%

plot(scrublet_adj)
```

![](doublet_detection_files/figure-html/scrublet%20manual%20threshold-1.png)

``` r

scrublet_adj_dt <- get_data(scrublet_adj)
scrublet_adj_dt[, Barcode := get_cell_names(sc_object)]
scrublet_adj_dt <- merge(scrublet_adj_dt, demuxlet_data, by = "Barcode")

doublet_metrics(
  predicted = scrublet_adj_dt$doublet,
  actual = scrublet_adj_dt$Call
)
#> $precision
#> [1] 0.6829566
#> 
#> $recall
#> [1] 0.7171264
#> 
#> $f1
#> [1] 0.6996245
```

We have made it worse here, but maybe there are cases where the scores
are not very good. Compared to the original version, Otsu’s method is
used to identify the perfect threshold, but it might fail in some
situations.

#### Running with grouping variables or within specific cell subsets

Scrublet (and the other methods) have the two parameters if you wish to
run the detection method within a specific set of cells (see the
`cells_to_use = NULL` parameter )and/or want to use a grouping variable
(like run the doublet detection within individual samples, see the
`group_by = NULL` parameter).

### Boosted doublet detection

The boosted approach trains a classifier to separate real cells from
synthetic doublets. This can work better when the Scrublet score
distribution is messy or unimodal. However, in this case we are running
25 iterations of the algorithm, which means 20x generation of synthetic
doublets. The initial PCA will be stored in the object for subsequent
reprojection - nonetheless, this makes this usually slower than the
other methods in this vignette.

``` r

boosted_doublets <- doublet_detection_boost_sc(
  object = sc_object,
  boost_params = params_boost()
)
```

Let’s check how boosted works?

``` r

boosted_dt <- get_data(boosted_doublets)
boosted_dt[, Barcode := get_cell_names(sc_object)]
boosted_dt <- merge(boosted_dt, demuxlet_data, by = "Barcode")

doublet_metrics(predicted = boosted_dt$doublet, actual = boosted_dt$Call)
#> $precision
#> [1] 0.6626324
#> 
#> $recall
#> [1] 0.5618987
#> 
#> $f1
#> [1] 0.6081222
```

Overall, it takes longer and has worse performance (in this data set).
If you wish to have a partially orthogonal method, it can be worthwhile
nonetheless. There is also an option to accelerate the algorithm via
fast clustering. This will run k-means clustering with
`sqrt(n_cells) * 4` centroids by default, run the kNN on these subset of
centroids, Louvain on the smaller kNN and project back the Louvain
membership for final scoring. However, this comes at cost of Recall, as
you can see below. Generally speaking, you need to see how many
centroids you need to represent your data properly. On smaller data
sets, you can get away with `n_cells / 10` and compress 10 cells into
one centroid… On larger data sets, you might want to compress further.

``` r

boosted_doublets_fast <- doublet_detection_boost_sc(
  object = sc_object,
  boost_params = params_boost(
    fast_cluster = TRUE,
    # in this case, we set the n_centroids to N_cells / 10
    # you want this to not be too low; also the data is so small, we just
    # use standard k-means over mini batch
    fast_cluster_params = list(n_centroids = 1500L, km_type = "standard")
  )
)
```

``` r

boosted_fast_dt <- get_data(boosted_doublets_fast)
boosted_fast_dt[, Barcode := get_cell_names(sc_object)]
boosted_fast_dt <- merge(boosted_fast_dt, demuxlet_data, by = "Barcode")

doublet_metrics(
  predicted = boosted_fast_dt$doublet,
  actual = boosted_fast_dt$Call
)
#> $precision
#> [1] 0.7209515
#> 
#> $recall
#> [1] 0.5054522
#> 
#> $f1
#> [1] 0.5942685
```

### scDblFinder

Lastly, `bixverse` also provides an implementation of the
[scDblFinder](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html)
from [Germain et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9204188/)
with some modifications: instead of xgboost, we will be using an
implementation of lightGBM using the histogram trick on the feature
space. This usually yields faster fitting compared to the original
xgboost approach. To run it, you can just do the following:

``` r

scdblfinder <- scdblfinder_sc(
  object = sc_object,
  scdblfinder_params = params_scdblfinder(
    expected_doublet_rate = 0.12
  ),
  # if you wish to get the features that were used to train the classifier
  return_features = FALSE
)
```

The scores are similar to the other methods.

``` r

scdblfinder_dt <- get_data(scdblfinder)

scdblfinder_dt[, Barcode := get_cell_names(sc_object)]
scdblfinder_dt <- merge(scdblfinder_dt, demuxlet_data, by = "Barcode")

doublet_metrics(
  predicted = scdblfinder_dt$predicted_doublets,
  actual = scdblfinder_dt$Call
)
#> $precision
#> [1] 0.592613
#> 
#> $recall
#> [1] 0.6895446
#> 
#> $f1
#> [1] 0.6374148
```

We can also extract the other scores from scDblFinder. We can for
example extract the weighted score, identify the threshold via Otsu’s
method and check the doublets like this.

``` r

weighted_scores <- get_scores(scdblfinder, score_type = "weighted")

weighted_scores_threshold <- find_threshold_otsu(x = weighted_scores)

weighted_doublets <- weighted_scores >= weighted_scores_threshold

table(
  lightgbm = scdblfinder_dt$predicted_doublets,
  weighted = weighted_doublets
)
#>         weighted
#> lightgbm FALSE  TRUE
#>    FALSE 12008   706
#>    TRUE     52  1762
```

Or alternatively, the cxds scores can also be extracted and used. These
scores were originally developed by [Bais and
Kostka](https://pmc.ncbi.nlm.nih.gov/articles/PMC7703774/) and are used
as an input into the classifier.

``` r

cxds_scores <- get_scores(scdblfinder, score_type = "cxds_scores")

cxds_scores_threshold <- find_threshold_otsu(x = cxds_scores)

cxds_doublets <- cxds_scores >= cxds_scores_threshold

table(
  lightgbm = scdblfinder_dt$predicted_doublets,
  cxds = cxds_doublets
)
#>         cxds
#> lightgbm FALSE  TRUE
#>    FALSE 12464   250
#>    TRUE    747  1067
```

Should you observe that this is very flat or you only have a few data
points (this can occur in low complexity experiments in which the
highest expressing genes are VERY similar across all the cells for
example in cell culture experiments), you can bump `cxds_genes` up. The
default here is 500L, but you might want to increase this pending your
experiment.

### Conclusion

All methods can help to identify in a data-driven way doublets in your
data. In this example, Scrublet performed better compared to the other
two methods; in a comparison paper scDblFinder is the best across
different benchmark data sets, see
[here](https://www.cell.com/cell-systems/fulltext/S2405-4712(20)30459-2).

## Clean up

``` r

unlink(tempdir_doublet, recursive = TRUE, force = TRUE)
```
