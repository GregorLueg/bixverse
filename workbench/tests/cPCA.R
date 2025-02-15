library(magrittr)
library(zeallot)

## Synthetic data ----

devtools::document()
devtools::load_all()

cPCA_data = create_synthetic_cPCA_data()

cPCA_data$target[1:5, 1:5]

## CoVar tests ----

ncol = 10000
nrow = 100

set.seed(123)
x = matrix(rnorm(ncol * nrow), nrow, ncol)

tictoc::tic()
y_1 = cov(x)
tictoc::toc()

tictoc::tic()
y_2 = rs_covariance(x)
tictoc::toc()

## Class tests ----

raw_data = t(cPCA_data$target)
background_mat = t(cPCA_data$background)

devtools::document()

sample_meta = data.table(
  sample_id = rownames(raw_data),
  grp = cPCA_data$target_labels,
  case_control = "case"
)

bulk_coexp_class = bulk_coexp(raw_data = raw_data, meta_data = sample_meta)

bulk_coexp_class = cPCA_preprocessing(bulk_coexp_class, background_mat = background_mat)

## Pre-processing ----

get_params(bulk_coexp_class, TRUE, TRUE)

target_mat <- S7::prop(bulk_coexp_class, "raw_data")
background_mat <- background_mat
verbose = TRUE
scale = FALSE

intersecting_features <- intersect(
  colnames(target_mat),
  colnames(background_mat)
)
if (verbose)
  message(sprintf(
    "A total of %i features/genes were identified",
    length(intersecting_features)
  )
  )

target_mat <- target_mat[, intersecting_features]
background_mat <- background_mat[, intersecting_features]

if (scale) {
  target_mat <- scale(target_mat, scale = scale)
  background_mat <- scale(background_mat, scale = scale)
}

target_covar = rs_covariance(target_mat)
background_covar = rs_covariance(background_mat)

dim(target_mat)

## cPCA functions ----

devtools::document()


alpha = 1
nPCs = 2

final_covar = target_covar - alpha * background_covar

final_covar[1:10, 1:10]

isSymmetric.matrix(final_covar)

tictoc::tic()
cPCA_results <- irlba::partial_eigen(final_covar, nPCs)
tictoc::toc()

cPCA_results$vectors

rextendr::document()

test_eigen(final_covar)

rs_contrastive_pca

c(loadings, factors) %<-% rs_contrastive_pca(
  target_covar = target_covar,
  background_covar = background_covar,
  target_mat = target_mat,
  alpha = .5,
  n_pcs = 5,
  return_loadings = FALSE
)

x = NULL

