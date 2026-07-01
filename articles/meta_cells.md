# Meta cells with bixverse

## Intro

Meta cells are small groups of transcriptionally similar cells
aggregated into a single representative profile. The motivation is
twofold: firstly, single cell counts are sparse and noisy, and many
downstream methods (correlation networks, GRN inference, archetype-based
modelling) behave poorly on raw single cell data. Aggregating into meta
cells reduces sparsity while preserving the heterogeneity that bulk
pseudo-replicates would average out. Secondly, they reduce the
computational burden substantially.

`bixverse` provides three meta cell algorithms with different underlying
ideas:

- **hdWGCNA-style meta cells** ([Morabito, et al.,
  2023](https://doi.org/10.1016/j.crmeth.2023.100498)) iterate over a
  kNN graph, picking seed cells and aggregating their neighbours, with a
  constraint on how many cells two meta cells are allowed to share.
  Simple, fast and reasonably effective. Works directly on a kNN graph.

- **SEACells** ([Persad, et al.,
  2023](https://doi.org/10.1038/s41587-023-01716-9)) uses kernel
  archetypal analysis. The algorithm finds a set of archetypes (the
  SEACells) such that every cell can be expressed as a convex
  combination of nearby archetypes. Tends to produce purer aggregations
  but is the most computationally expensive of the three.

- **SuperCells** ([Bilous, et al.,
  2022](https://doi.org/10.1186/s12859-022-04861-1)) applies the
  walktrap community detection algorithm to the kNN graph and treats
  each community as a meta cell. Granularity is controlled indirectly
  through a graining factor (the average number of cells per meta cell).

One implementation detail worth flagging upfront. In `bixverse`, the
`SingleCells` class keeps counts on disk in DuckDB and binary Rust-based
files and streams them in as needed. The `MetaCells` class does not.
After aggregation the count matrix is small enough to hold entirely in
memory as a sparse matrix, and all downstream operations (HVG, PCA,
neighbours, clustering, embeddings) run against this in-memory
representation (most of the code still makes usage of Rust’s incredible
performance nonetheless). Practically that means the meta cell pipeline
is substantially faster than the equivalent operations on the parent
`SingleCells` object for two reasons: data is already in memory and more
importantly, it is simply less data.

``` r

library(bixverse)
library(ggplot2)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
library(patchwork)
library(magrittr)
```

## Preparing the single cell data

We use the same CD34 cells from the [SEACells
vignette](https://github.com/dpeerlab/SEACells/blob/main/notebooks/SEACell_computation.ipynb)
to set up the parent `SingleCells` object: load (QC not needed - it’s
already filtered data), HVG selection, PCA, kNN graph. We can use the
provided cell type labels to check purity.

``` r

cd34_path <- download_cd34_data()

tempdir_cd34 <- tempdir()

sc_object <- SingleCells(dir_data = tempdir_cd34)

sc_object <- load_h5ad(object = sc_object, h5_path = cd34_path)
#>  Using light streaming for the CSR to CSC conversion.
#> Loading observations data from h5ad into the DuckDB.
#> Loading variables data from h5ad into the DuckDB.
```

Let’s run HVG detection, PCA and kNN (+ sNN) generation.

``` r

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = 2000L,
  .verbose = FALSE
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 30L,
  sparse_svd = FALSE
)
#> Using dense SVD solving on scaled data on 2000 HVG.

sc_object <- find_neighbours_sc(
  object = sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(ann_dist = "euclidean", knn_method = "kmknn")
  )
)
#> 
#> Generating sNN graph (full: TRUE).
#> Transforming sNN data to igraph.

# we need the kNN object for the diffusion maps
knn_object <- get_knn_obj(x = sc_object)
```

## Generating meta cells

All three generators take a `SingleCells` object, read counts from the
binary files as needed, and return a `MetaCells` object. They share a
common output structure: an obs table with one row per meta cell
(recording which original cells went into it), a var table mirrored from
the parent object, and an in-memory sparse count matrix.

### hdWGCNA-style (bootstrapped metacells)

The hdWGCNA-style algorithm reuses the kNN graph already on the object.
The two parameters that matter most are `target_no_metacells` and
`max_shared`, which controls how many cells two meta cells may have in
common. The original paper uses `1000L` meta cells, but we want to also
compare the different metrics downstream, so, we will set it to `250L`.
We also recalculate the kNN graph with larger k, to include more cells
per given meta cell.

``` r

hdwgcna <- generate_bt_meta_cells_sc(
  object = sc_object,
  sc_meta_cell_params = params_sc_bt_metacells(
    target_no_metacells = 250L,
    knn = list(k = 25L)
  ),
  regenerate_knn = TRUE, # regenerate kNN graph with 25 neighbours
  .verbose = TRUE
)

hdwgcna
#> Single cell experiment (Meta Cells).
#>   Meta cell method: meta_cells_hdwgcna
#>   No meta cells: 250
#>   No genes: 12464
#>   No original cells: 6881
#>   No unassigned cells: 2977
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

Importantly, not all cells HAVE to be assigned to a meta cell and cells
can occur in several meta cells at once. We can appreciate here that a
number of cells remains unassigned and each MetaCell does contain self +
neighbours, i.e., 26 (`self + k = 25L`) original cells. Let’s check how
often the cells occur… ?

``` r

no_duplicated_cells <- table(unlist(hdwgcna[[]]$original_cell_idx))

hist(
  no_duplicated_cells,
  xlab = "No of times a cell occurs",
  main = "No of times a cell is part of a meta cell",
  breaks = 10L
)
```

![](meta_cells_files/figure-html/hdwgcna%20-%20cells%20per%20meta%20cell-1.png)

As we can appreciate, some of the cells are indeed shared across meta
cells. Let’s add the diffusion coordinates to this one for later
analysis

``` r

hdwgcna <- calc_diffusion_coordinates(object = hdwgcna, knn_data = knn_object)
```

### SEACells

SEACells operates on the PCA embedding (or any batch corrected
embedding) and the kNN graph and runs an iterative kernel archetypal
analysis. `n_sea_cells` sets the number of archetypes;
`min_iter`/`max_iter` bound the optimisation; `convergence_epsilon` sets
the early stopping threshold relative to the initial RSS. The algorithm
was aggressively optimised to be faster. With `pruning = TRUE` small
values are set to zero, reducing calculations… The rule of thumb by the
authors is one SEACell per 75 cells, i.e., 90 SEACells for this
experiment. We will set this higher here to `250L` to make this a bit
more comparable with the other methods…

``` r

seacells <- generate_seacells_sc(
  object = sc_object,
  seacell_params = params_sc_seacells(
    n_sea_cells = 250L,
    min_iter = 10L,
    convergence_epsilon = 0.001,
    pruning = TRUE
  ),
  .verbose = TRUE
)

seacells
#> Single cell experiment (Meta Cells).
#>   Meta cell method: seacell
#>   No meta cells: 250
#>   No genes: 12464
#>   No original cells: 6881
#>   No unassigned cells: 0
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

We can appreciate that each SEACell can contain a varying degree of
original cells contributing to this specific metacell.

``` r

no_cells_seacells <- purrr::map_dbl(seacells[[]]$original_cell_idx, length)

hist(
  no_cells_seacells,
  xlab = "No cells per meta cell",
  main = "Cells per SEACells meta cell",
  breaks = 25L
)
```

![](meta_cells_files/figure-html/seacells%20-%20cells%20per%20meta%20cell-1.png)

We will also add the diffusion coordinates to this one

``` r

seacells <- calc_diffusion_coordinates(object = seacells, knn_data = knn_object)
```

#### SEACells on large data

SEACells is quite infamous for being very slow and difficult to run on
larger data sets. bixverse offers some tricks here which make the
generation of SEACells feasible on large data sets:

1.  The initialisation is the part of the original algorithm that scales
    the worst. The greedy column subset selection needs a full `K^2`
    column for every candidate cell, which is `O(N2)` and simply will
    not finish on millions of cells. `bixverse` therefore switches
    strategy by size: above `greedy_threshold` it drops the greedy
    top-up and samples the archetypes at random, which is effectively
    free. When a better-than-random start is wanted but the greedy pass
    is still out of reach, setting `n_landmarks` enables a Nyström
    route: a small set of density-weighted landmarks is chosen, the
    diffusion operator is built and eigendecomposed on those landmarks
    alone, and the multiscale embedding is then carried to all cells via
    a Nyström extension before the usual max-min waypoint sampling. The
    expensive eigendecomposition stays at landmark scale `L x L` instead
    of `N x N`.
2.  The optimisation loop itself runs in bounded memory. `K^2` is never
    formed; every `K^2 * X` term is evaluated as `K (K x X)` so memory
    stays bounded by the number of non-zeros in `K` rather `N^2`. The
    Frank-Wolfe gradients are computed one column at a time and in
    parallel, so the full `k x N` gradient never has to exist at once,
    and pruning keeps A and B sparse across iterations.
3.  The RSS for the convergence check is computed differently above
    20,000 cells. Rather than materialising the `N x N` reconstruction
    `K x B x A`, the squared residual is expanded with the trace
    identity and evaluated through cyclic reordering, so every
    intermediate is at most `N x k` or `k x k`. The value is the same
    Frobenius norm, just without the dense reconstruction. Together
    these are what made it possible to run on a million cells locally,
    which the original implementation cannot do.

If you want to run SEACells on larger data sets, these are the knobs you
have to make this run in a reasonable time on a powerful laptop.

### SuperCells

SuperCells runs walktrap on the kNN graph. The number of meta cells is
set indirectly through `graining_factor`: with a factor of 30 over ~6800
cells you should get roughly ~230 meta cells.

``` r

supercells <- generate_supercells_sc(
  object = sc_object,
  sc_supercell_params = params_sc_supercell(
    graining_factor = 30
  ),
  .verbose = TRUE
)

supercells
#> Single cell experiment (Meta Cells).
#>   Meta cell method: supercells
#>   No meta cells: 230
#>   No genes: 12464
#>   No original cells: 6881
#>   No unassigned cells: 0
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

Let’s plot the number of cells per SuperCell.

``` r

no_cells_supercells <- purrr::map_dbl(supercells[[]]$original_cell_idx, length)

hist(
  no_cells_supercells,
  xlab = "No cells per meta cell",
  main = "Cells per Supercell",
  breaks = 25L
)
```

![](meta_cells_files/figure-html/supercells%20-%20cells%20per%20meta%20cell-1.png)

Similar to SEACells, we get a gradient here.

``` r

supercells <- calc_diffusion_coordinates(
  object = supercells,
  knn_data = knn_object
)
```

#### SuperCells on large data

SuperCell is already far cheaper than SEACells on large data: there is
no archetypal analysis, just Walktrap community detection on the kNN
graph. The one part that does not scale naively is the random-walk
representation. `compute_walk_probabilities()` within the Rust code
gives every cell a sparse vector of landing probabilities for a
`walk_length`-step walk; on a well-connected graph these vectors spread
to touch a large fraction of cells after only a few steps, so the store
drifts towards dense and memory grows as `O(n × support)`.

`max_support` caps this. With `max_support = k` each initial walk vector
keeps only its `k` largest entries (the dropped tail holds negligible
walk mass), bounding the store at roughly `k × n` regardless of how far
the walks spread. This makes the run an approximation: the Ward merge
criterion is driven by distances between walk vectors, and truncating
them can shift the merge order slightly, so you may get marginally
different communities than the exact run. With `max_support = NULL` (the
default) the walks are kept exact. Beyond saving memory, a smaller `k`
also speeds up every distance and merge operation, since both are linear
in the support.

``` r

supercells_large <- generate_supercells_sc(
  object = sc_object,
  sc_supercell_params = params_sc_supercell(
    graining_factor = 30,
    max_support = 256L # would bound the walk vectors to 256
  ),
  .verbose = TRUE
)
```

## Metrics

### Purity

A simple sanity check is to ask, for each meta cell, what fraction of
its constituent cells share the same underlying label. We use the Leiden
clusters from the parent `SingleCells` object as a proxy for cell type
identity. The labels themselves are imperfect, so absolute numbers
should be read with that caveat, but relative differences between
methods are reasonably indicative.

``` r

cell_labels <- unlist(sc_object[["celltype"]]) %>% as.character()

hdwgcna <- calc_meta_cell_purity(hdwgcna, original_cell_type = cell_labels)
seacells <- calc_meta_cell_purity(seacells, original_cell_type = cell_labels)
supercells <- calc_meta_cell_purity(
  supercells,
  original_cell_type = cell_labels
)

purity_dt <- rbind(
  data.table(method = "hdWGCNA", purity = hdwgcna[[]]$mc_purity),
  data.table(method = "SEACells", purity = seacells[[]]$mc_purity),
  data.table(method = "SuperCells", purity = supercells[[]]$mc_purity)
)

purity_dt[, .(mean = mean(purity), median = median(purity)), by = method]
#>        method      mean    median
#>        <char>     <num>     <num>
#> 1:    hdWGCNA 0.8947692 1.0000000
#> 2:   SEACells 0.9025032 1.0000000
#> 3: SuperCells 0.8979310 0.9806216
```

``` r

ggplot(purity_dt, aes(x = method, y = purity)) +
  geom_violin(aes(fill = method), alpha = 0.6) +
  geom_boxplot(width = 0.05, outlier.size = 0.5) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Method", y = "Meta cell purity (Leiden)") +
  ylim(0, 1)
```

![](meta_cells_files/figure-html/purity%20plot-1.png)

In terms of cell type purity, the methods are basically the same.
Potentially more interesting is their behaviour in the manifold, see
below…

### Manifold regions

We called the
[`calc_diffusion_coordinates()`](https://gregorlueg.github.io/bixverse/reference/calc_diffusion_coordinates.md)
with the kNN graph on the data. This will tell us which parts of the
manifold are being sampled here, i.e., for example high, medium or low
density regions (based on the diffusion map and distance to the 150k-th
neighbour within that). This gives an idea how many rare cell states the
meta cells capture.

``` r

hdwgcna_regions <- as.data.table(
  table(hdwgcna[[]]$density_region) / nrow(hdwgcna[[]])
)[, method := "hdWGCNA"]
seacells_regions <- as.data.table(
  table(seacells[[]]$density_region) / nrow(seacells[[]])
)[, method := "SEACells"]
supercells_regions <- as.data.table(
  table(supercells[[]]$density_region) / nrow(supercells[[]])
)[, method := "SuperCells"]

region_dt <- rbind(
  hdwgcna_regions,
  seacells_regions,
  supercells_regions
)[, V1 := factor(V1, levels = c("high", "mid", "low"))]
setnames(region_dt, old = c("V1", "N"), new = c("region", "proportion"))
```

``` r

ggplot(region_dt, aes(x = method, y = proportion, fill = region)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Method", y = "Proportion of meta cells", fill = "Density region") +
  scale_fill_manual(
    values = setNames(
      c("#2c2d54", "#969bc7", "#6f9954"),
      c("high", "mid", "low")
    )
  )
```

![](meta_cells_files/figure-html/regions%20plot-1.png)

This is one of the most important differences between the methods.
hdWGCNA will disproportionally sample regions of the manifolds with
medium density and basically proportion-based. SEACells on the other
hand generates meta cells that capture sparse (more heterogenous)
regions of the manifold. SuperCells sits between the two methods.

### Compactness and separation

We can also assess other metrics on the manifold representation.
Compactness (how close are the cells within a given metacell to the
centroid) and separation (how far is the closest centroid). The first
gives an (unbiased) indication of purity; the latter of diversity
captured in the manifold. We can do this trivially via:

``` r

hdwgcna <- calc_manifold_metrics(hdwgcna)
seacells <- calc_manifold_metrics(seacells)
supercells <- calc_manifold_metrics(supercells)
```

Let’s generate a plotting data.table

``` r

metrics_dt <- rbind(
  data.table(
    method = "hdWGCNA",
    separation = hdwgcna[[]]$separation,
    compactness = hdwgcna[[]]$compactness,
    region = hdwgcna[[]]$density_region
  ),
  data.table(
    method = "SEACells",
    separation = seacells[[]]$separation,
    compactness = seacells[[]]$compactness,
    region = seacells[[]]$density_region
  ),
  data.table(
    method = "SuperCells",
    separation = supercells[[]]$separation,
    compactness = supercells[[]]$compactness,
    region = supercells[[]]$density_region
  )
)[, region := factor(region, levels = c("high", "mid", "low"))]

metrics_dt[,
  .(
    mean_separation = mean(separation),
    median_separation = median(separation),
    mean_compactness = mean(compactness),
    median_compactness = median(compactness)
  ),
  .(method)
]
#>        method mean_separation median_separation mean_compactness
#>        <char>           <num>             <num>            <num>
#> 1:    hdWGCNA       0.1255286        0.08895643      0.007053792
#> 2:   SEACells       0.1802643        0.15110685      0.011203591
#> 3: SuperCells       0.2016617        0.15913163      0.006335426
#>    median_compactness
#>                 <num>
#> 1:        0.004067266
#> 2:        0.004984010
#> 3:        0.003216350
```

Some differences are visible here. hdWGCNA has the worst separation and
best compactness (not unexpected given the bootstrapping method only
sampling direct neighbours). SEACells has the worst compactness across
the three methods (not unsurprising as the meta cells here capture more
sparse regions of the manifold) and middling separation; SuperCells
(version 2.0, the default setting here) is in the middle.

``` r

per_region_stats <- metrics_dt[,
  .(
    mean_separation = mean(separation),
    median_separation = median(separation),
    mean_compactness = mean(compactness),
    median_compactness = median(compactness)
  ),
  .(method, region)
]
setorder(per_region_stats, method, region)

per_region_stats[]
#>        method region mean_separation median_separation mean_compactness
#>        <char> <fctr>           <num>             <num>            <num>
#> 1:   SEACells   high      0.04633167        0.05101920     0.0007073069
#> 2:   SEACells    mid      0.13877805        0.13594331     0.0044319796
#> 3:   SEACells    low      0.24090206        0.21325481     0.0186658600
#> 4: SuperCells   high      0.04950358        0.05186578     0.0006019381
#> 5: SuperCells    mid      0.14154094        0.13623749     0.0031448643
#> 6: SuperCells    low      0.30667216        0.25426194     0.0113716719
#> 7:    hdWGCNA   high      0.03334880        0.02947950     0.0009821354
#> 8:    hdWGCNA    mid      0.10448298        0.09349236     0.0055345098
#> 9:    hdWGCNA    low      0.26961609        0.19062942     0.0168331578
#>    median_compactness
#>                 <num>
#> 1:       0.0006614107
#> 2:       0.0033688664
#> 3:       0.0076578967
#> 4:       0.0005725418
#> 5:       0.0024108913
#> 6:       0.0055384417
#> 7:       0.0009728607
#> 8:       0.0042686835
#> 9:       0.0093960930
```

We can also look at this per density region in the manifold where we can
observe that compactness is smallest in the high density regions,
whereas separation is the highest in the low density regions. And below
as plots:

#### Compactness plots

``` r

ggplot(metrics_dt, aes(x = method, y = compactness)) +
  geom_violin(aes(fill = method), alpha = 0.6) +
  geom_boxplot(width = 0.05, outlier.size = 0.5) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Method", y = "Compactness")
```

![](meta_cells_files/figure-html/compactness%20plot-1.png)

``` r

ggplot(metrics_dt, aes(x = method, y = compactness)) +
  geom_boxplot(aes(fill = region), outlier.size = 0.5) +
  theme_bw() +
  labs(x = "Method", y = "Compactness") +
  scale_fill_manual(
    values = setNames(
      c("#2c2d54", "#969bc7", "#6f9954"),
      c("high", "mid", "low")
    )
  )
```

![](meta_cells_files/figure-html/compactness%20per%20region-1.png)

#### Separation plots

``` r

ggplot(metrics_dt, aes(x = method, y = separation)) +
  geom_violin(aes(fill = method), alpha = 0.6) +
  geom_boxplot(width = 0.05, outlier.size = 0.5) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Method", y = "Separation")
```

![](meta_cells_files/figure-html/separation%20plot-1.png)

``` r

ggplot(metrics_dt, aes(x = method, y = separation)) +
  geom_boxplot(aes(fill = region), outlier.size = 0.5) +
  theme_bw() +
  labs(x = "Method", y = "Compactness") +
  scale_fill_manual(
    values = setNames(
      c("#2c2d54", "#969bc7", "#6f9954"),
      c("high", "mid", "low")
    )
  )
```

![](meta_cells_files/figure-html/separation%20per%20region-1.png)

Based on this, we can conclude the following:

- Purity is VERY similar across the methods (for this data set).
- At matched K (~250 metacells), SuperCells achieves the lowest median
  compactness, slightly ahead of hdWGCNA and SEACells. Separation is
  comparable between SuperCells and SEACells, with hdWGCNA noticeably
  lower.
- The methods differ markedly in where they place metacells. SEACells
  assigns 50% of metacells to low-density regions versus 24% for hdWGCNA
  and 37% for SuperCells, consistent with its kernel archetypal
  formulation seeking out underrepresented states.
- Method choice should reflect the downstream question. For tight
  aggregation in dense regions (e.g., bulk-like pseudo-replicates of
  common cell types), SuperCells is efficient and effective. For
  preserving rare populations or transitioning states (e.g.,
  differentiation trajectories), SEACells’ rare-state bias is desirable,
  with the trade-off of more heterogeneous metacells in those regions.
  hdWGCNA is fast and simple but assignments overlap more in embedding
  space.

## Working with the MetaCells class

This is where the in-memory representation pays off. HVG selection, PCA,
neighbour graphs and embeddings on a few hundred meta cells take seconds
and do not touch the Rust binary files on disk. We will just take
forward SEACells here, but the methods below work across all of
`MetaCells`.

### HVGs and PCA

The same HVG and PCA method dispatches you know from single cells work
here…

``` r

seacells <- find_hvg_sc(
  object = seacells,
  hvg_no = 2000L,
  .verbose = FALSE
)

seacells <- calculate_pca_sc(
  object = seacells,
  no_pcs = 30L
)
```

### Neighbours, clustering and UMAP

And also the same interfaces for neighbours, (Leiden) clustering and
UMAP can be used here.

``` r

seacells <- find_neighbours_sc(
  object = seacells,
  neighbours_params = params_sc_neighbours(
    knn = list(k = 10L, knn_method = "exhaustive")
  )
)
#> 
#> Generating sNN graph (full: TRUE).
#> Transforming sNN data to igraph.

seacells <- find_clusters_sc(seacells, res = 1.0, name = "leiden_clusters")

seacells <- umap_sc(seacells, k = 10L, knn_method = "exhaustive")
#> Running UMAP.
#> Using n_epochs = 500 (dataset <10k samples or adam_parallel optimiser)
#> Using provided kNN graph.
```

``` r

umap_dt <- as.data.table(
  get_embedding(seacells, "umap"),
  keep.rownames = "meta_cell_id"
)[, leiden_clusters := seacells[["leiden_clusters"]]]

ggplot(umap_dt, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(colour = as.factor(leiden_clusters)), size = 1.5) +
  theme_bw() +
  labs(colour = "MetaCell cluster")
```

![](meta_cells_files/figure-html/mc%20umap%20plot-1.png)

The plot is much sparser than the equivalent UMAP on the original cells,
which is the point: each dot is a denoised aggregate of ~15-20 cells,
and the structure that survives aggregation is the structure that’s
robust to sparsity.

### Co-expression module detection on meta cells

Meta cells compress information and reduce sparsity. That makes them
ideal for co-expression module detection methods and we can use a SCENIC
implementation within bixverse [akin to the single cell
version](https://gregorlueg.github.io/bixverse/articles/analysis_single_cell.html#scenic)
also on MetaCells.

Let’s run this quickly…

``` r

tf_dt <- data.table::fread(
  "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt",
  header = FALSE,
  col.names = "tf"
)

scenic_res <- scenic_grn_sc(
  object = seacells,
  tf_ids = tf_dt$tf,
  scenic_params = params_scenic(
    learner_type = "randomforest",
    gene_batch_size = 64L,
    # due to the small data set size, we should set 'min_samples_leaf' lower
    # than for a massive single cell data set
    learner_params = list(min_samples_leaf = 10L)
  ),
  .verbose = TRUE
)
#> No target genes supplied, running gene filter...
#> SCENIC gene filter: 12464 / 12464 genes pass.
#> Warning in `method(scenic_grn_sc, bixverse::MetaCells)`(object = <object>, :
#> 610 TF identifier(s) not found in the object and dropped.
#> SCENIC: 12464 target genes, 1282 TFs, 250 cells
```

``` r

scenic_res <- identify_tf_to_genes(
  scenic_res,
  n_sd = 2,
  .verbose = TRUE
)
#> Extracting TF to gene associations via per-gene threshold (mean + 2.0 * SD).

scenic_res <- tf_to_genes_correlations(
  x = scenic_res,
  object = seacells,
  cor_filter = 0.01,
  .verbose = TRUE
)
#> Calculating the pairwise correlations between the TFs and genes
#> Removing TF <> gene pairs with cors <= 0.010
#> Removing self loops (TF controlling its own expression

tf_gene_dt <- get_tf_to_gene(scenic_res)

# we will not do any filtering
tf_to_gene_ls <- split(tf_gene_dt$gene, tf_gene_dt$tf)
```

In a proper situation we would filter down the TF to gene associations
via motifs, please refer to the details in this
[vignette](https://gregorlueg.github.io/bixverse/articles/analysis_single_cell.html#scenic).
We will skip this step for now and just run AUCell to show how it works.

``` r

auc_res <- aucell_sc(
  object = seacells,
  gs_list = tf_to_gene_ls
)

umap_dt[, `:=`(IRF1 = auc_res[, "IRF1"], TCF4 = auc_res[, "TCF4"])]

p1 <- ggplot(umap_dt, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(fill = IRF1), size = 2.5, shape = 21) +
  theme_bw() +
  labs(fill = "IRF1 AUC") +
  scale_fill_viridis_c()

p2 <- ggplot(umap_dt, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(fill = TCF4), size = 2.5, shape = 21) +
  theme_bw() +
  labs(fill = "TCF4 AUC") +
  scale_fill_viridis_c()

p1 + p2
```

![](meta_cells_files/figure-html/run%20AUCell%20on%20meta%20cells-1.png)

## Pseudo-bulking

Pseudo-bulks are conceptually adjacent but solve a different problem.
Where meta cells are dense aggregates designed to feed into
co-expression and archetype methods, pseudo-bulks are coarse summaries
(typically one bulk per sample-by-cluster combination) used as input to
bulk-style DGE tools such as limma-voom, edgeR or DESeq2 to get around
pseudo-replication problems and p-value inflation in cell-based DGEs,
see [Zimmerman et
al.](https://www.nature.com/articles/s41467-021-21038-1)

`get_pseudobulked_sc` takes a named list of cell IDs and returns either
a dense matrix or a sparse `dgRMatrix`. With `assay = "raw"` it sums raw
counts (what bulk DGE tools expect); with `assay = "norm"` it averages
normalised counts.

``` r

cl_dt <- sc_object[[c("leiden", "cell_id")]]
cell_list <- split(cl_dt$cell_id, cl_dt$leiden)

pb_counts <- get_pseudobulked_sc(
  object = sc_object,
  cell_list = cell_list,
  return_format = "sparse",
  assay = "raw",
  .verbose = FALSE
)

dim(pb_counts)
#> [1]    14 12464
```

The returned matrix has one row per group and one column per gene. From
here the standard bulk DGE pipeline applies; `bixverse` does not wrap
the DGE call itself.

## Clean up

``` r

unlink(tempdir_cd34, recursive = TRUE, force = TRUE)
```
