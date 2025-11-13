# Contrastive PCA

## Contrastive PCA

`bixverse` provides an implementation of the contrastive PCA algorithm
from [Abid, et al.](https://www.nature.com/articles/s41467-018-04608-8).
This vignette will show you how to run the algorithm with some synthetic
data.

``` r
library(bixverse)
library(data.table)
```

### Data

Let’s explore first the synthetic data (same as from the original
authors).

``` r
cpca_test_data <- synthetic_c_pca_data()

plot(cpca_test_data)
```

![](cpca_files/figure-html/synthetic%20data-1.png)

Synthetic data visualization for Contrastive PCA analysis

The pre-dominant signal is in the first third of the features (noise)
basically making it quite difficult to distinguish the 4 groups in the
target data. The background data reproduces the same pre-dominant
signal, but lacks the more subtle signals. This makes it ideal to remove
the background signal and enhance the more salient features.

### Explore the alpha parameter space

Let’s prepare the data for further downstream analysis. As these are to
some degree a type of co-expression module analysis, we can generate the
corresponding object. These objects are based on
[S7](https://github.com/RConsortium/S7).

``` r
raw_data <- t(cpca_test_data$target)
background_mat <- t(cpca_test_data$background)

sample_meta <- data.table(
  sample_id = rownames(raw_data),
  grp = cpca_test_data$target_labels
)

cpca_obj <- bulk_coexp(raw_data = raw_data, meta_data = sample_meta)

# This function allows you to do scaling, HVG selection, etc.
# We won't be doing any of that
cpca_obj <- preprocess_bulk_coexp(cpca_obj)
#> A total of 30 genes will be included.
```

Now we need to prepare the object for the application of contrastive
PCA. To do so, we will provide the background data.

``` r
cpca_obj <- contrastive_pca_processing(
  cpca_obj,
  background_mat = background_mat
)
#> A total of 30 features/genes were identified
```

The object is now ready to explore the alpha parameter. Alpha = 0 is
equivalent to standard PCA. `bixverse` provides a function to quickly
plot out the impact of different alphas.

``` r
c_pca_plot_alphas(
  cpca_obj,
  label_column = "grp",
  n_alphas = 10L,
  max_alpha = 100
)
#> Found the grp in the meta data. Adding labels to the graph.
```

![](cpca_files/figure-html/plot%20different%20alphas-1.png)

Impact of the alpha parameter on the ability to distinguish the groups

In the range of alpha ~1.5 to alpha 10 we can appreciate that the four
different groups are now visible, whereas at alpha 0 everything is just
a blob. Once you have chosen an alpha parameter, you can run:

### Run contrastive PCA and extract the data

``` r
# run contrastive pca
cpca_obj <- contrastive_pca(
  cpca_obj,
  alpha = 2.5,
  no_pcs = 10L
)

# extract the loadings
c_pca_loadings <- get_c_pca_loadings(cpca_obj)

# extract the factors
c_pca_factors <- get_c_pca_factors(cpca_obj)
```

The data can then be used for subsequent further analysis.
